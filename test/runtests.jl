using OceanScalingTests
using Test
using MPI

@testset "OceanScalingTests.jl" begin
   println("Hello world")
end

function run_doubledrake_simulation(ranks; experiment = :DoubleDrake)
   simulation = OceanScalingTests.scaling_test_simulation(1/5, ranks, 10minutes, stop_time; 
                                                          child_arch = CPU(), Nz = 10, experiment)

   run!(simulation)

   return true
end

@testset "DoubleDrake runs" begin
   MPI.Init()
   
   comm   = MPI.COMM_WORLD
   rank   = MPI.Comm_rank(comm)
   Nranks = MPI.Comm_size(comm)
   
   Ry = 1
   Rx = 1
   
   ranks = (Rx, Ry, 1)
   
   @test run_doubledrake_simulation(ranks)
end
   
