using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using MPI
using JLD2

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

ranks      = (Nranks, 1, 1)
resolution = 6

bathymetry = jldopen("/home/ssilvest/development/OceanScalingTests.jl/data/bathymetry_six.jld2")["bathymetry"]

Δt = 5minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, rank, bathymetry, Δt, stop_iteration)

MPI.finalize()