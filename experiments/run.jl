using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using NVTX
using MPI
using JLD2

using Oceananigans.Models.HydrostaticFreeSurfaceModels

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

ranks       = (Rx, Ry, 1)
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "DoubleDrake"))

restart = get(ENV, "RESTART", "")
restart_iteration = "28800"


Δt = 10minutes * (3 / resolution)
stop_time = 100days

Δt = 35

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_time experiment restart
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_time; experiment, restart)

if !isnothing(simulation)
    OceanScalingTests.set_outputs!(simulation, Val(experiment), overwrite_existing = true)
    
    if isempy(restart)
        run!(simulation)
    else
        run!(simulation, pickup = "restart/RealisticOcean_checkpoint_$(rank)_iteration$(restart_iteration).jld2")
    end

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

# MPI.Finalize()
