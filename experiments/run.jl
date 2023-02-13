using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using MPI
using JLD2

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

ranks       = (Rx, Ry, 1)
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "Quiescent"))
use_buffers = parse(Bool, get(ENV, "USEBUFFERS", "1"))

Δt = 10minutes * (3 / resolution)
stop_iteration = 4000

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_iteration experiment use_buffers
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_iteration; experiment, use_buffers)

OceanScalingTests.set_outputs!(simulation, Val(experiment))

run!(simulation)

# MPI.Finalize()
