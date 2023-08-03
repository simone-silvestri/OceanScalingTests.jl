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
output_dir     = get(ENV, "OUTPUTDIR", "./")

final_year  = parse(Int, get(ENV, "FINALYEAR",  "0"))
final_month = parse(Int, get(ENV, "FINALMONTH", "12"))

max_Δt = precision(45 * 48 * 1.5 / resolution)
min_Δt = 5 # precision(45 * 48 * 0.9 / resolution)
stop_time = 7300days

if rank == 0
    @info "Scaling test" ranks resolution max_Δt min_Δt stop_time experiment profile with_fluxes with_restoring restart 
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, (min_Δt, max_Δt), stop_time; Nz, experiment, restart,
						       profile, with_fluxes, with_restoring, loadbalance, precision)

if !isnothing(simulation)
    @info "type of dt :" typeof(simulation.Δt)
    OceanScalingTests.set_outputs!(simulation, Val(experiment); output_dir, overwrite_existing = false, checkpoint_time = 10days)

    pickup_file = isempty(restart) ? false :  output_dir * "RealisticOcean_checkpoint_$(rank)_iteration$(restart).jld2"
    run!(simulation, pickup = pickup_file)

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

MPI.Finalize()
