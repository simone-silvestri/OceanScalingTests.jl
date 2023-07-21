using Preferences
const iscray = parse(Bool, load_preference(Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"), "iscray", "false"))
@debug "Preloading GTL library" iscray
if iscray
    import Libdl
    Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
end
using MPI
MPI.Init()

using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime, SpecifiedTimes
using NVTX
using JLD2

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

@show Rx, Ry

ranks       = (Rx, Ry, 1)

# Enviromental variables
resolution     = parse(Int,  get(ENV, "RESOLUTION", "3"))
experiment     = Symbol(     get(ENV, "EXPERIMENT", "DoubleDrake"))
with_fluxes    = parse(Bool, get(ENV, "WITHFLUXES", "0"))
with_restoring = parse(Bool, get(ENV, "WITHRESTORING", "0"))
profile        = parse(Bool, get(ENV, "PROFILE", "1"))
restart        =             get(ENV, "RESTART", "")
Nz             = parse(Int,  get(ENV, "NZ", "100"))
loadbalance    = parse(Bool, get(ENV, "LOADBALANCE", "0"))
precision      = eval(Symbol(get(ENV, "PRECISION", "Float64")))

final_year  = parse(Int, get(ENV, "FINALYEAR",  "0"))
final_month = parse(Int, get(ENV, "FINALMONTH", "12"))

Δt = 20 # precision(45 * 48 / resolution)
stop_time = 7300days

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_time experiment profile with_fluxes with_restoring restart 
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_time; Nz, experiment, restart,
						       profile, with_fluxes, with_restoring, loadbalance, precision)

increase_Δt1(sim) = sim.Δt = 30
increase_Δt2(sim) = sim.Δt = 40
increase_Δt3(sim) = sim.Δt = precision(45 * 48 / resolution)

if !isnothing(simulation)
    @info "type of dt :" typeof(simulation.Δt)
    OceanScalingTests.set_outputs!(simulation, Val(experiment); overwrite_existing = true, checkpoint_time = 10days)

    simulation.callbacks[:increase_Δt1] = Callback(increase_Δt1, SpecifiedTimes([60days]))
    simulation.callbacks[:increase_Δt2] = Callback(increase_Δt2, SpecifiedTimes([120days]))
    simulation.callbacks[:increase_Δt3] = Callback(increase_Δt3, SpecifiedTimes([150days]))

    pickup_file = isempty(restart) ? false :  "./RealisticOcean_checkpoint_$(rank)_iteration$(restart).jld2"
    run!(simulation, pickup = pickup_file)

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

MPI.Finalize()
