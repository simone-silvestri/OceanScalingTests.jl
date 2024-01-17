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
                             CFL = 0.75)
    wave_speed = sqrt(g * grid.Lz)
    local_Δ    = 1 / sqrt(1 / minimum_xspacing(grid)^2 + 1 / minimum_yspacing(grid)^2)
    global_Δ   = MPI.Allreduce(local_Δ, min, grid.architecture.communicator)

    return max(Int(ceil(2 * Δt / (CFL / wave_speed * global_Δ))), 10)
end

substeps(free_surface) = length(free_surface.settings.substepping.averaging_weights)
# At the moment the simulation does not run for TEOS10 with the DoubleDrake initialization (srt(negative number))
equation_of_state(::Val{E},            precision) where E = TEOS10EquationOfState(precision)
equation_of_state(::Val{:DoubleDrake}, precision)         = LinearEquationOfState(precision)

experiment_depth(E) = E == :RealisticOcean ? 5244.5 : 3000

using Oceananigans.Advection: CrossAndSelfUpwinding, EnergyConserving

best_momentum_advection(grid, precision) = VectorInvariant(vorticity_scheme = WENO(precision; order = 9),
							    vertical_scheme = Centered(),
							  divergence_scheme = WENO(precision))

simple_momentum_advection(grid, precision) = VectorInvariant()

function scaling_test_simulation(resolution, ranks, Δt, stop_time;
                                 child_arch = GPU(),
                                 experiment = :Quiescent, 
                                 Depth = experiment_depth(experiment),
                                 latitude = (-75, 75),
                                 restart = "",
                                 z_faces_function = exponential_z_faces,
                                 κ_skew = 0,
                                 κ_symmetric = 0,
                                 Nz = 100,
                                 profile = false,
                                 with_fluxes = true,
                                 with_restoring = true,
                                 loadbalance = true,
                                 slope_limiter = FluxTapering(1e-2),
                                 precision = Float64,
                                 boundary_layer_parameterization = RiBasedVerticalDiffusivity(precision)
                                 )

    arch = Distributed(child_arch; partition = Partition(ranks...))

    min_Δt, max_Δt = Δt isa Number ? (Δt, Δt) : Δt

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
    
    horizontal_closure = if κ_skew > 0 || κ_symmetric > 0
        IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric, slope_limiter)
    else
        nothing
    end

    closure  = (vertical_diffusivity, boundary_layer_parameterization, horizontal_closure)

    tracer_advection   = Oceananigans.Advection.ThreeDimensionalTracerAdvection(;
                                x = WENO(precision; order = 7),
                                y = WENO(precision; order = 7),
                                z = Centered(precision))
    
    momentum_advection = best_momentum_advection(grid, precision)

    free_surface = SplitExplicitFreeSurface(precision; substeps = barotropic_substeps(max_Δt, grid))

    @info "running with $(substeps(free_surface)) barotropic substeps"

    buoyancy = SeawaterBuoyancy(precision; equation_of_state=equation_of_state(Val(experiment), precision))
    
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
                                          momentum_advection, 
                                          tracer_advection,
                                          coriolis,
                                          buoyancy,
                                          tracers,
                                          boundary_conditions,
                                          closure)

    @info "model allocated"

    #####
    ##### Initial condition:
    #####

    # If we are profiling launch only 100 time steps and mark each one with NVTX
    if profile
        profiled_time_steps!(model, max_Δt, resolution)
        return nothing
    end
   
    initialize_model!(model, Val(experiment); restart)
    @info "model initialized"

    #####
    ##### Simulation
    #####

    simulation = Simulation(model; Δt=min_Δt, stop_time)
 
    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u
        v = sim.model.velocities.v
        w = sim.model.velocities.w
        η = sim.model.free_surface.η
	T = sim.model.tracers.T
	S = sim.model.tracers.S

        rk = sim.model.grid.architecture.local_rank

	    @info @sprintf("R: %02d, T: % 12s, it: %d, Δt: %.2f, vels: %.2e %.2e %.2e, η: %.2e, trac: %.2e %.2e, wt : %s", 
                        rk, prettytime(sim.model.clock.time), sim.model.clock.iteration, sim.Δt, 
                        maximum(abs, u),  maximum(abs, v), maximum(abs, w), maximum(abs, η), maximum(abs, T), maximum(abs, S),
                        prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    wizard = TimeStepWizard(cfl=0.35; max_change=1.1, max_Δt, min_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    
    if experiment == :RealisticOcean 
        if with_fluxes 
            simulation.callbacks[:update_fluxes] = Callback(update_fluxes!, TimeInterval(5days))
        end
        if with_restoring 
            repeat_year = parse(Bool, get(ENV, "REPEATYEAR", "true"))
            update_time = repeat_year ? 1days : 15days
            simulation.callbacks[:update_restoring] = Callback(update_restoring!, TimeInterval(update_time))
        end
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
