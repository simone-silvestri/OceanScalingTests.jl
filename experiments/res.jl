using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using MPI
using JLD2

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = Nranks
Rx = 1

ranks       = (Rx, Ry, 1)
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "Quiescent"))
use_buffers = parse(Bool, get(ENV, "USEBUFFERS", "0"))

Δt = 15minutes * (3 / resolution)
stop_iteration = 5000

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_iteration experiment
end

OceanScalingTests.run_scaling_test!(resolution, ranks, Δt, stop_iteration; experiment, use_buffers)

MPI.Finalize()