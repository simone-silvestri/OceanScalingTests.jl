using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.Utils 
using Oceananigans.Units
using Oceananigans.Architectures: device
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using SeawaterPolynomials: TEOS10EquationOfState
using JLD2

# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid; 
                             g = Oceananigans.BuoyancyModels.g_Earth, 
                             CFL = 0.7)
    wave_speed = sqrt(g * grid.Lz)
    local_Δ    = 1 / sqrt(1 / minimum_xspacing(grid)^2 + 1 / minimum_yspacing(grid)^2)
    global_Δ   = MPI.Allreduce(local_Δ, min, grid.architecture.communicator)

    return max(Int(ceil(2 * Δt / (CFL / wave_speed * global_Δ))), 10)
end

substeps(free_surface) = length(free_surface.settings.substepping.averaging_weights)
# At the moment the simulation does not run for TEOS10 with the DoubleDrake initialization (srt(negative number))
equation_of_state(::Val{E},            precision) where E = TEOS10EquationOfState(precision)
equation_of_state(::Val{:DoubleDrake}, precision)         = LinearEquationOfState(precision)

using Oceananigans.Advection: CrossAndSelfUpwinding, EnergyConservingScheme

previous_momentum_advection(grid, precision) = VectorInvariant(vorticity_scheme = WENO(precision),
                                                                vertical_scheme = WENO(grid),
                                                             ke_gradient_scheme = EnergyConservingScheme(precision),
                                                                      upwinding = CrossAndSelfUpwinding()) 

experiment_depth(E) = E == :RealisticOcean ? 5244.5 : 3000

best_momentum_advection(grid, precision) = VectorInvariant(vorticity_scheme = WENO(precision; order = 9),
                                                            vertical_scheme = WENO(grid))

function scaling_test_simulation(resolution, ranks, Δt, stop_time;
                                 child_arch = GPU(),
                                 experiment = :Quiescent, 
                                 Depth = experiment_depth(experiment),
                                 latitude = (-75, 75),
                                 restart = "",
                                 z_faces_function = exponential_z_faces,
                                 Nz = 100,
                                 profile = false,
                                 with_fluxes = true,
                                 with_restoring = true,
                                 loadbalance = true,
                                 precision = Float64,
                                 boundary_layer_parameterization = RiBasedVerticalDiffusivity(precision),
                                 diffuse_initially = true)

    topo = (Periodic, Bounded, Bounded)
    arch = DistributedArch(child_arch; topology = topo, ranks)

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
    
    tracer_advection   = WENO(grid)
    momentum_advection = best_momentum_advection(grid, precision)

    free_surface = SplitExplicitFreeSurface(precision; substeps = barotropic_substeps(Δt, grid))

    @info "running with $(substeps(free_surface)) barotropic substeps"

    buoyancy = SeawaterBuoyancy(precision; equation_of_state=equation_of_state(Val(experiment), precision))
    closure  = (vertical_diffusivity, boundary_layer_parameterization)
    
    if diffuse_initially
        closure_init = (closure..., HorizontalScalarDiffusivity(κ = 100))
    else
        closure_init = closure
    end

    coriolis = HydrostaticSphericalCoriolis(precision)

    #####
    ##### Boundary conditions
    #####

    boundary_conditions = set_boundary_conditions(Val(experiment), grid; with_fluxes, with_restoring)

    #####
    ##### Model setup
    #####

    tracers = boundary_layer_parameterization isa CATKEVerticalDiffusivity ?
              (:T, :S, :e) : (:T, :S)

    @info "allocating model"
    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection = nothing, tracer_advection = nothing,
                                          coriolis = nothing,
                                          buoyancy,
                                          tracers,
                                          closure = closure_init)

    @info "model allocated"

    #####
    ##### Initial condition:
    #####

    # If we are profiling launch only 100 time steps and mark each one with NVTX
    if profile
        profiled_time_steps!(model, Δt, resolution)
        return nothing
    end
   
    initialize_model!(model, Val(experiment); restart)
    @info "model initialized"

    #####
    ##### Simulation setup
    #####

    stop_iteration = 1000    
    simulation = Simulation(model; Δt, stop_iteration)
    run!(simulation)

    T, S = model.tracers

    @info "reallocating model"
    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection, tracer_advection,
                                          coriolis,
                                          buoyancy,
                                          tracers,
                                          boundary_conditions,
                                          closure)

    @info "model reallocated"

    set!(model, T = T, S = S)    
    @info "model reinitialized"

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
            simulation.callbacks[:update_fluxes] = Callback(update_fluxes!, TimeInterval(5days))
        end
        if with_restoring 
            simulation.callbacks[:update_restoring] = Callback(update_restoring!, TimeInterval(15days))
        end
        ## simulation.callbacks[:garbage_collect] = Callback((sim) -> GC.gc(), IterationInterval(50))
    end

    return simulation
end

function profiled_time_steps!(model, Δt, resolution; gc_steps = 10, profiled_steps = 30)
    # initial time steps
    for step in 1:10
        time_step!(model, Δt)
    end
    
    nranks = MPI.Comm_size(MPI.COMM_WORLD)
    rank   = MPI.Comm_rank(MPI.COMM_WORLD)

    if rank == 0
      @info "start profiling"
    end
    elapsed_time = Float64[0] 
    # Perform profiling
    for step in 1:gc_steps
        for nogc in 1:profiled_steps
            elapsed_time[1] += @elapsed begin
	            NVTX.@range "one time step" begin
                    time_step!(model, Δt)
                end
            end
        end
    end

    elapsed_time[1] = elapsed_time[1] / nranks / gc_steps / profiled_steps
    MPI.Allreduce!(elapsed_time, +, MPI.COMM_WORLD)

    if rank == 0
	file = "time_res$(resolution)_NZ$(model.grid.Nz)_ranks$(nranks)_prec$(eltype(model.grid)).jld2"
        while isfile(file)
	        file = "new_" * file
        end
        jldsave(file, elapsed_time=elapsed_time[1])
    end

    MPI.Barrier(MPI.COMM_WORLD)

    return nothing
end
