using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Test
using MPI
using Pkg

@testset "OceanScalingTests.jl" begin
   println("Hello world")
end

function coarse_simulation_runs(experiment, ranks)
   simulation = OceanScalingTests.scaling_test_simulation(1/5, ranks, 10minutes, 365days; 
                                                          child_arch = CPU(), Nz = 10, experiment)

   simulation.stop_iteration = 10
   run!(simulation)

   return true
end

@testset "Quiescent and DoubleDrake run" begin
   MPI.Init()
   
   comm   = MPI.COMM_WORLD
   rank   = MPI.Comm_rank(comm)
   Nranks = MPI.Comm_size(comm)
   
   Ry = 1
   Rx = 1
   
   ranks = (Rx, Ry, 1)
   
   @test coarse_simulation_runs(:Quiescent,   ranks)
   @test coarse_simulation_runs(:DoubleDrake, ranks)
end

@testset "Generate RealisticOcean fluxes" begin
   
   ENV["FINALYEAR"] = "1996"
   include("../fluxes/generate_fluxes.jl")

   generate_fluxes(1; arch = CPU())
   fluxes = filter(x -> x[end-3:end] == "jld2", readdir("../fluxes"))

   @test length(fluxes) == 
end