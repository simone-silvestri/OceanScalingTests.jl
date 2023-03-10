import Pkg
Pkg.activate("/nobackup/users/ssilvest/OceanScalingTests.jl/")
using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using NVTX
using MPI
using JLD2

using Oceananigans.Models.HydrostaticFreeSurfaceModels

MPI.Init()

include("/nobackup/users/ssilvest/OceanScalingTests.jl/precompile_0/precompile_Oceananigans.Models.HydrostaticFreeSurfaceModels.jl")
_precompile_()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

ranks       = (Rx, Ry, 1)
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "DoubleDrake"))
use_buffers = parse(Bool, get(ENV, "USEBUFFERS", "1"))

Δt = 10minutes * (3 / resolution)
stop_iteration = 30000

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_iteration experiment use_buffers
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_iteration; experiment, use_buffers)

if !isnothing(simulation)
    OceanScalingTests.set_outputs!(simulation, Val(experiment))
    run!(simulation)

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

# MPI.Finalize()
