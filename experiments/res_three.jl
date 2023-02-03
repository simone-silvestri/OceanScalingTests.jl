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

ranks      = (Rx, Ry, 1)
resolution = 3

bathymetry = jldopen("/home/ssilvest/development/OceanScalingTests.jl/data/bathymetry_three.jld2")["bathymetry"]

Δt = 5minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, Δt, stop_iteration; bathymetry, latitude = (-75, 75))

MPI.finalize()