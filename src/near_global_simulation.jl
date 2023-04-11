using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.Utils 
using Oceananigans.Units
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using SeawaterPolynomials: TEOS10EquationOfState
using CUDA

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
                                 boundary_layer_parameterization = RiBasedVerticalDiffusivity(),
                                 Nz = 100,
                                 profile = false,
                                 with_fluxes = true)

    child_arch = GPU()

    topo = (Periodic, Bounded, Bounded)
    arch = DistributedArch(child_arch; topology = topo, ranks, use_buffers)

    Lφ = latitude[2] - latitude[1]

    # grid size
    Nx = Int(360 * resolution)
    Ny = Int(Lφ * resolution)

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
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(realistic_bathymetry(underlying_grid))) :
           experiment == :DoubleDrake ?
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(double_drake_bathymetry)) :
           underlying_grid

    #####
    ##### Physics setup and numerical methods
    #####

    νz = 5e-4
    κz = 3e-5        

    vertical_diffusivity = VerticalScalarDiffusivity(ν=νz, κ=κz)
    horizont_diffusivity = HorizontalScalarDiffusivity(ν=ad_hoc_viscosity, discrete_form=true, 
                                                       loc = (Center, Center, Center),
                                                       parameters = (sᵐᵃˣ = 2.5, νᶜ = 1000.0))
    
    tracer_advection   = WENO(underlying_grid)
    momentum_advection = VectorInvariant(vorticity_scheme  = WENO(), 
                                         divergence_scheme = WENO(), 
                                         vertical_scheme   = WENO(underlying_grid)) 

    free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(Δt, grid, g_Earth))

    @info "running with $(free_surface.settings.substeps) barotropic substeps"

    buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
    closure  = (vertical_diffusivity, boundary_layer_parameterization, horizont_diffusivity)
    
    coriolis = HydrostaticSphericalCoriolis()

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

    initialize_model!(model, Val(experiment); restart)
    @info "model initialized"

    # If we are profiling launch only 100 time steps and mark each one with NVTX
    if profile
        profiled_time_step!(model, Δt)
        return nothing
    end

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

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))
    
    if experiment == :RealisticOcean 
        if with_fluxes
            simulation.callbacks[:update_fluxes]   = Callback(update_fluxes, TimeInterval(5days))
        end
        simulation.callbacks[:garbage_collect] = Callback((sim) -> GC.gc(), IterationInterval(50))
    end

    return simulation
end

const metrics = [
    "sm__cycles_elapsed.avg",
    "sm__cycles_elapsed.avg.per_second",

    "dram__bytes.sum",
    "lts__t_bytes.sum",
    "l1tex__t_bytes.sum",

    "sm__sass_thread_inst_executed_op_fadd_pred_on.sum",
    "sm__sass_thread_inst_executed_op_fmul_pred_on.sum",
    "sm__sass_thread_inst_executed_op_ffma_pred_on.sum",

    "sm__sass_thread_inst_executed_op_dadd_pred_on.sum",
    "sm__sass_thread_inst_executed_op_dmul_pred_on.sum",
    "sm__sass_thread_inst_executed_op_dfma_pred_on.sum",

    "sm__sass_thread_inst_executed_op_hadd_pred_on.sum",
    "sm__sass_thread_inst_executed_op_hmul_pred_on.sum",
    "sm__sass_thread_inst_executed_op_hfma_pred_on.sum",
]

function process(metrics)
    time = metrics["sm__cycles_elapsed.avg"] / metrics["sm__cycles_elapsed.avg.per_second"]
    D_FLOP = 2*metrics["sm__sass_thread_inst_executed_op_dfma_pred_on.sum"] + 
               metrics["sm__sass_thread_inst_executed_op_dmul_pred_on.sum"] +
               metrics["sm__sass_thread_inst_executed_op_dadd_pred_on.sum"]
    F_FLOP = 2*metrics["sm__sass_thread_inst_executed_op_ffma_pred_on.sum"] + 
               metrics["sm__sass_thread_inst_executed_op_fmul_pred_on.sum"] +
               metrics["sm__sass_thread_inst_executed_op_fadd_pred_on.sum"]
    H_FLOP = 2*metrics["sm__sass_thread_inst_executed_op_hfma_pred_on.sum"] + 
               metrics["sm__sass_thread_inst_executed_op_hmul_pred_on.sum"] +
               metrics["sm__sass_thread_inst_executed_op_hadd_pred_on.sum"]

    AI_D_DRAM = D_FLOP / metrics["dram__bytes.sum"]
    AI_F_DRAM = F_FLOP / metrics["dram__bytes.sum"]
    AI_H_DRAM = H_FLOP / metrics["dram__bytes.sum"]

    AI_D_L2 = D_FLOP / metrics["lts__t_bytes.sum"]
    AI_F_L2 = F_FLOP / metrics["lts__t_bytes.sum"]
    AI_H_L2 = H_FLOP / metrics["lts__t_bytes.sum"]

    AI_D_L1 = D_FLOP / metrics["l1tex__t_bytes.sum"]
    AI_F_L1 = F_FLOP / metrics["l1tex__t_bytes.sum"]
    AI_H_L1 = H_FLOP / metrics["l1tex__t_bytes.sum"]

    D_FLOPs = D_FLOP/time
    H_FLOPs = H_FLOP/time
    F_FLOPs = F_FLOP/time

    @info "Kernel performance" time D_FLOP F_FLOP H_FLOP D_FLOPs F_FLOPs H_FLOPs
    @info "Arithmetic intensity (DRAM)" AI_D_DRAM AI_F_DRAM AI_H_DRAM
    @info "Arithmetic intensity (L2)" AI_D_L2 AI_F_L2 AI_H_L2
    @info "Arithmetic intensity (L1)" AI_D_L1 AI_F_L1 AI_H_L1
end

using NVPerfWorks


function profiled_time_step!(model, Δt; gc_steps = 10, profiled_steps = 10)
    # initial time steps
    for step in 1:10
        time_step!(model, Δt)
    end

    measure = NVPerfWorks.StatefulMeasure(metrics)
  
    # Perform profiling
    for step in 1:gc_steps
        for nogc in 1:profiled_steps
            NVPerfWorks.start!(measure)
            stats = @timed time_step!(model, Δt)
            measured_metrics = NVPerfWorks.stop!(measure)
            @info "Duration" time=stats.time gctime=stats.gctime bytes=stats.bytes
            if measured_metrics === nothing
                continue
            end
           
            @info "Fraction active"  fraction = measure.sessions/measure.passes
            process(measured_metrics)
        end
        GC.gc()
    end
end
