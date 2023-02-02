using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using MPI
using JLD2

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Rx = Ry = Int(sqrt(Nranks))

ranks      = (Rx, Ry, 1)
resolution = 8

bathymetry = jldopen("/home/ssilvest/development/OceanScalingTests.jl/data/bathymetry_eight.jld2")["bathymetry"]

Δt = 3.5minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, rank, bathymetry, Δt, stop_iteration)

MPI.finalize()