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

# Enviromental variables
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "DoubleDrake"))
with_fluxes = parse(Bool, get(ENV, "WITHFLUXES", "1"))
profile     = parse(Bool, get(ENV, "PROFILE", "1"))
restart     = get(ENV, "RESTART", "")


Δt = 10minutes * (3 / resolution)
stop_time = 100days

Δt = 40
Nz = 100

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_time experiment profile with_fluxes restart 
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_time; Nz, experiment, restart, profile, with_fluxes)

if !isnothing(simulation)
    OceanScalingTests.set_outputs!(simulation, Val(experiment); overwrite_existing = true, checkpoint_time = 10days)
    
    if isempty(restart)
        run!(simulation)
    else
        pickup_file =  "restart/RealisticOcean_checkpoint_$(rank)_iteration$(restart).jld2"
        @info "restarting from $(pickup_file)"
        run!(simulation, pickup = pickup_file)
    end

    @info "simulation took $(prettytime(simulation.run_wall_time))"
end

# MPI.Finalize()
