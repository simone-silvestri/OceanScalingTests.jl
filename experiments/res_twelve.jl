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

ranks      = (Rx, Ry, 1)
resolution = 12

bathymetry = nothing # jldopen("/home/ssilvest/development/OceanScalingTests.jl/data/bathymetry_twelve.jld2")["bathymetry"]

Δt = 2.5minutes
stop_iteration = 5000

OceanScalingTests.run_scaling_test!(resolution, ranks, rank, Δt, stop_iteration; bathymetry)

MPI.finalize()
