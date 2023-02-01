using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using MPI

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

ranks      = (Nranks, 1, 1)
resolution = 12

bathymetry = jldopen("../data/bathymetry_twelve.jld2")["bathymetry"]

Δt = 2.5minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, bathymetry, Δt, stop_iteration)

MPI.finalize()