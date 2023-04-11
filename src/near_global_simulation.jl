using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Utils 
using Oceananigans.Units
using Oceananigans.Architectures: device
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using SeawaterPolynomials: TEOS10EquationOfState
using CUDA: synchronize

smoothing_convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz=10.0, convective_νz=10.0,
                                                                          background_κz=1.0,  background_νz=1.0)

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)
    global_Δ   = MPI.Allreduce(local_Δ, min, grid.architecture.communicator)

    return max(Int(ceil(2 * Δt / (CFL / wave_speed * global_Δ))), 10)
end

@inline function ad_hoc_viscosity(i, j, k, grid, clock, fields, p) 
    speed = spᶜᶜᶜ(i, j, k, grid, fields)
    return ifelse(speed > p.sᵐᵃˣ, p.νᶜ, zero(grid))
end

experiment_depth(exp) = exp == :RealisticOcean ? 5244.5 : 3kilometers

function scaling_test_simulation(resolution, ranks, Δt, stop_time;
                                 experiment = :Quiescent, 
                                 Depth = experiment_depth(experiment),
                                 latitude = (-75, 75),
                                 use_buffers = true,
                                 restart = "",
                                 z_faces_function = exponential_z_faces,
                                 Nz = 100,
                                 profile = false,
                                 with_fluxes = true,
                                 loadbalance = true,
				 precision = Float64,
				 boundary_layer_parameterization = RiBasedVerticalDiffusivity(precision))

    child_arch = GPU()

    topo = (Periodic, Bounded, Bounded)
    arch = DistributedArch(child_arch; topology = topo, ranks, use_buffers)

    Lφ = latitude[2] - latitude[1]

    # grid size
    Nx = Int(360 * resolution) 
    Ny = Int(Lφ * resolution)

    z_faces = z_faces_function(Nz, Depth)

    grid = load_balanced_grid(arch, precision, (Nx, Ny, Nz), latitude, z_faces, resolution, Val(loadbalance), Val(experiment))

    #####
    ##### Physics setup and numerical methods
    #####

    νz = 5e-4
    κz = 3e-5        

    vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), precision; ν=νz, κ=κz)
    horizont_diffusivity = HorizontalScalarDiffusivity(precision; ν=ad_hoc_viscosity, discrete_form=true, 
                                                       loc = (Center, Center, Center),
						       parameters = (sᵐᵃˣ = precision(2.5), νᶜ = precision(1000.0)))
    
    tracer_advection   = WENO(grid.underlying_grid)
    momentum_advection = VectorInvariant(vorticity_scheme  = WENO(precision), 
                                         divergence_scheme = WENO(precision), 
					 vertical_scheme   = WENO(grid.underlying_grid))

    free_surface = SplitExplicitFreeSurface(; gravitational_acceleration = precision(g_Earth),
					      substeps = barotropic_substeps(Δt, grid, g_Earth))

    @info "running with $(free_surface.settings.substeps) barotropic substeps"

    buoyancy = SeawaterBuoyancy(precision; equation_of_state=TEOS10EquationOfState(precision))
    closure  = (vertical_diffusivity, boundary_layer_parameterization) #, horizont_diffusivity)
    
    coriolis = HydrostaticSphericalCoriolis(precision)

    #####
    ##### Boundary conditions
    #####

    boundary_conditions = set_boundary_conditions(Val(experiment), grid; with_fluxes)

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
 
    # If we are profiling launch only 100 time steps and mark each one with NVTX
    if profile
        profiled_time_steps!(model, Δt)
        return nothing
    end

    initialize_model!(model, Val(experiment); restart)
    @info "model initialized"

    #####
    ##### Simulation setup
    #####

    simulation = Simulation(model; Δt, stop_time)
 
    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u
        v = sim.model.velocities.v
        w = sim.model.velocities.w
        η = sim.model.free_surface.η

        rk = sim.model.grid.architecture.local_rank

	    @info @sprintf("Rank: %02d, Time: % 12s, it: %d, Δt: %.2f, vels: %.2e ms⁻¹ %.2e ms⁻¹ %.2e ms⁻¹, max(|η|): %.2e m, wt : %s", 
                        rk, prettytime(sim.model.clock.time), sim.model.clock.iteration, sim.Δt, 
                        maximum(abs, u),  maximum(abs, v), maximum(abs, w), maximum(abs, η),
                        prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    
    if experiment == :RealisticOcean 
        if with_fluxes
            simulation.callbacks[:update_fluxes] = Callback(update_fluxes, TimeInterval(5days))
        end
        simulation.callbacks[:garbage_collect] = Callback((sim) -> Threads.@spawn GC.gc(), IterationInterval(50))
    end

    return simulation
end

function profiled_time_steps!(model, Δt; gc_steps = 50, profiled_steps = 30)
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
	#synchronize()
	#MPI.Barrier(MPI.COMM_WORLD)
        GC.gc()
    end
end
