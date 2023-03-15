using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Utils 
using Oceananigans.Units
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)
    global_Δ   = MPI.Allreduce(local_Δ, min, grid.architecture.communicator)

    return Int(ceil(2 * Δt / (CFL / wave_speed * global_Δ)))
end

experiment_depth(exp) = exp == :RealisticOcean ? 6kilometers : 3kilometers

function scaling_test_simulation(resolution, ranks, Δt, stop_iteration;
                                 Depth = experiment_depth(experiment),
                                 experiment = :Quiescent, 
                                 latitude = (-75, 75),
                                 use_buffers = false,
                                 z_faces_function = exponential_z_faces,
                                 boundary_layer_parameterization = RiBasedVerticalDiffusivity())

    child_arch = GPU()

    topo = (Periodic, Bounded, Bounded)
    arch = DistributedArch(child_arch; topology = topo, ranks, use_buffers)

    Lφ = latitude[2] - latitude[1]

    # grid size
    Nx = Int(360 * resolution)
    Ny = Int(Lφ * resolution)
    Nz = 120

    z_faces = z_faces_function(Nz, Depth)

    # A spherical domain
    @show underlying_grid = LatitudeLongitudeGrid(arch,
                                                  size = (Nx, Ny, Nz),
                                                  longitude = (-180, 180),
                                                  latitude = latitude,
                                                  halo = (5, 5, 5),
                                                  z = z_faces,
                                                  precompute_metrics = true)

    grid = experiment == :RealisticOcean ? 
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(realistic_bathymetry(grid))) :
           experiment == :DoubleDrake ?
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(double_drake_bathymetry)) :
           underlying_grid

    #####
    ##### Physics setup and numerical methods
    #####

    νz = 5e-4
    κz = 3e-5

    vertical_diffusivity  = VerticalScalarDiffusivity(ν=νz, κ=κz)
        
    tracer_advection   = WENO(underlying_grid)
    momentum_advection = VectorInvariant(vorticity_scheme  = WENO(), 
                                         divergence_scheme = WENO(), 
                                         vertical_scheme   = WENO(underlying_grid)) 

    free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(Δt, grid, g_Earth))

    @info "running with $(free_surface.settings.substeps) barotropic substeps"

    buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())
    closure  = (vertical_diffusivity, boundary_layer_parameterization)
    coriolis = HydrostaticSphericalCoriolis()

    #####
    ##### Boundary conditions
    #####

    boundary_conditions = set_boundary_conditions(Val(experiment), grid.Ly)

    #####
    ##### Model setup
    #####

    tracers = boundary_layer_parameterization isa CATKEVerticalDiffusivity ?
              (:T, :S, :e) : (:T, :S)

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection, tracer_advection,
                                          coriolis,
                                          buoyancy,
                                          tracers,
                                          boundary_conditions,
                                          closure)

    #####
    ##### Initial condition:
    #####

    initialize_model!(model, Val(experiment))
    @info "model initialized"

    @show model.velocities.u.boundary_conditions

    profile = parse(Bool, get(ENV, "PROFILE", "1"))
    
    # If we are profiling launch only 100 time steps and mark each one with NVTX
    if profile
        profiled_time_step!(model, Δt)
        return nothing
    end

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

    if experiment == :RealisticOcean
        simulation.callbacks[:update_fluxes] = Callback(update_fluxes, TimeInterval(9days))
    end

    return simulation
end

function profiled_time_step!(model, Δt; gc_steps = 10, profiled_steps = 10)
    # initial time steps
    for step in 1:10
        time_step!(model, Δt)
    end
  
    # Perform profiling
    for step in 1:gc_steps
        for nogc in 1:profiled_steps
            NVTX.@range "one time step" begin
                time_step!(model, Δt)
            end
        end
        GC.gc()
    end
end
