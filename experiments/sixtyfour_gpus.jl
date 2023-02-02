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
resolution = 24

bathymetry = jldopen("../data/bathymetry_twentyfour.jld2")["bathymetry"]

Δt = 1minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, rank, bathymetry, Δt, stop_iteration)

MPI.finalize()