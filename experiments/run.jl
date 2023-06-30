using Preferences
const iscray = parse(Bool, load_preference(Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"), "iscray", "false"))
@debug "Preloading GTL library" iscray
if iscray
    import Libdl
    Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
end
using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using NVTX
using MPI
using JLD2

using Oceananigans.Models.HydrostaticFreeSurfaceModels

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

@show Rx, Ry

ranks       = (Rx, Ry, 1)

# Enviromental variables
resolution  = parse(Int,  get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(     get(ENV, "EXPERIMENT", "DoubleDrake"))
with_fluxes = parse(Bool, get(ENV, "WITHFLUXES", "0"))
profile     = parse(Bool, get(ENV, "PROFILE", "1"))
restart     =             get(ENV, "RESTART", "")
Nz          = parse(Int,  get(ENV, "NZ", "100"))
loadbalance = parse(Bool, get(ENV, "LOADBALANCE", "0"))
precision   = eval(Symbol(get(ENV, "PRECISION", "Float64")))

Δt = precision(45 * 48 / resolution)
stop_time = 1000

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_time experiment profile with_fluxes restart 
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_time; Nz, experiment, restart,
						       profile, with_fluxes, loadbalance, precision)

if !isnothing(simulation)
    @info "type of dt :" typeof(simulation.Δt)
    OceanScalingTests.set_outputs!(simulation, Val(experiment); overwrite_existing = true, checkpoint_time = 10days)
    
    pickup = isempty(restart) ? false :  "./RealisticOcean_checkpoint_$(rank)_iteration$(restart).jld2"
    run!(simulation, pickup = pickup_file)

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

MPI.Finalize()
