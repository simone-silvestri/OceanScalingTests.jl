using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Utils 
using Oceananigans.Distributed: partition_global_array

function run_scaling_test!(resolution, ranks, rank, bathymetry, Δt, stop_iteration;
                           no_ibg = false)

    child_arch = GPU()

    topo = (Periodic, Bounded, Bounded)
    arch = MultiArch(child_arch; topology = topo, ranks)

    latitude = (-75, 75)

    # 0.25 degree resolution
    Nx = 360 * resolution
    Ny = 150 * resolution
    Nz = 250

    N = (Nx, Ny, Nz)

    z_faces = linear_z_faces(Nz, Depth)

    rx, ry, _ = arch.local_index

    # A spherical domain
    @show underlying_grid = LatitudeLongitudeGrid(arch,
                                                  size = (Nx, Ny, Nz),
                                                  longitude = (-180, 180),
                                                  latitude = latitude,
                                                  halo = (5, 5, 5),
                                                  z = z_faces,
                                                  precompute_metrics = true)

    nx, ny, nz = size(underlying_grid)
    bathymetry = bathymetry[1 + nx * (rx - 1) : rx * nx, 1 + ny * (ry - 1) : ry * ny]

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

    #####
    ##### Physics and model setup
    #####

    νz = 5e-3
    κz = 1e-4

    convective_adjustment  = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.2, convective_νz = 0.2)
    vertical_diffusivity   = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=νz, κ=κz)
        
    tracer_advection   = WENO()
    momentum_advection = VectorInvariant(vorticity_scheme  = WENO(), 
                                         divergence_scheme = WENO(), 
                                         vertical_scheme   = WENO()) 

    #####

    CFL            = 0.8
    wave_speed     = sqrt(g_Earth * grid.Lz)
    Δg             = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)
    @show substeps = Int(ceil(2 * Δt / (CFL / wave_speed * Δg)))

    free_surface = SplitExplicitFreeSurface(; substeps)
    buoyancy     = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())
    
    closure      = (vertical_diffusivity, convective_adjustment)
    coriolis     = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection, tracer_advection,
                                          coriolis,
                                          buoyancy,
                                          tracers = (:T, :S),
                                          closure,
                                          calculate_only_active_cells_tendencies = no_ibg)

    #####
    ##### Initial condition:
    #####

    @info "model initialized"

    #####
    ##### Simulation setup
    #####

    simulation = Simulation(model; Δt, stop_iteration)

    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u

        @show @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wall time: %s", 
                        prettytime(sim.model.clock.time),
                        sim.model.clock.iteration, maximum(abs, u), 
                        prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    # Let's goo!
    @show "Running with Δt = $(prettytime(simulation.Δt))"

    run!(simulation)

    @show """
        Simulation took $(prettytime(simulation.run_wall_time))
        Free surface: $(typeof(model.free_surface).name.wrapper)
        Time step: $(prettytime(Δt))
    """
end
