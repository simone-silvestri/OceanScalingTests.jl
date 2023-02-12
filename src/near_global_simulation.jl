using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Utils 
using Oceananigans.Units
using Oceananigans.Distributed: partition_global_array

function scaling_test_simulation(resolution, ranks, Δt, stop_iteration;
                                 Depth = 3kilometers,
                                 experiment = :Quiescent, 
                                 latitude = (-80, 80),
                                 use_buffers = false,
                                 z_faces_function = exponential_z_faces)

    child_arch = GPU()

    topo = (Periodic, Bounded, Bounded)
    arch = DistributedArch(child_arch; topology = topo, ranks, use_buffers)

    Lφ = latitude[2] - latitude[1]

    # 0.25 degree resolution
    Nx = Int(360 * resolution)
    Ny = Int(Lφ * resolution)
    Nz = 150

    @show Ny / arch.ranks[2]

    z_faces = z_faces_function(Nz, Depth)

    # A spherical domain
    @show grid = LatitudeLongitudeGrid(arch,
                                       size = (Nx, Ny, Nz),
                                       longitude = (-180, 180),
                                       latitude = latitude,
                                       halo = (5, 5, 5),
                                       z = z_faces,
                                       precompute_metrics = true)


    if experiment == :DoubleDrake
        grid = ImmersedBoundaryGrid(grid, GridFittedBottom(double_drake_bathymetry))

	@show arch.local_rank, grid.immersed_boundary
    end

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
    CFL            = 0.6
    wave_speed     = sqrt(g_Earth * grid.Lz)
    Δgr            = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)

    Δg = MPI.Allreduce(Δgr, min, arch.communicator)

    @show substeps = Int(ceil(2 * Δt / (CFL / wave_speed * Δg)))

    free_surface = SplitExplicitFreeSurface(; substeps)
    buoyancy     = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())
    
    closure      = (vertical_diffusivity, convective_adjustment)
    coriolis     = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())

    boundary_conditions = set_boundary_conditions(Val(experiment))

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection, tracer_advection,
                                          coriolis,
                                          buoyancy,
                                          tracers = (:T, :S),
                                          boundary_conditions,
                                          closure)

    #####
    ##### Initial condition:
    #####

    initialize_model!(model, Val(experiment))
    @info "model initialized"

    #####
    ##### Simulation setup
    #####

    simulation = Simulation(model; Δt, stop_iteration)

    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u
        v = sim.model.velocities.v
        w = sim.model.velocities.w
        η = sim.model.free_surface.η

	@info @sprintf("Time: % 12s, iteration: %d, max(|u|, |v|, |w|): %.2e ms⁻¹ %.2e ms⁻¹ %.2e ms⁻¹, max(|η|): %.2e m, wall time: %s", 
                        prettytime(sim.model.clock.time),
                        sim.model.clock.iteration, maximum(abs, u),  maximum(abs, v), maximum(abs, w), maximum(abs, η),
                        prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    # Let's goo!
    @show "Running with Δt = $(prettytime(simulation.Δt))"

    return simulation
end
