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

ranks       = (Rx, Ry, 1)
resolution  = parse(Int, get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(get(ENV, "EXPERIMENT", "Quiescent"))
use_buffers = parse(Bool, get(ENV, "USEBUFFERS", "0"))

Δt = 10minutes * (3 / resolution)
stop_iteration = 100

if rank == 0
    @info "Scaling test" ranks resolution Δt stop_iteration experiment use_buffers
end

simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, Δt, stop_iteration; experiment, use_buffers)

model   = simulation.model

outpus  = Dict()
indices = (:, :, model.grid.Nz) 

outputs[:u] = Field(model.velocities.u; indices)
outputs[:v] = Field(model.velocities.v; indices)
outputs[:w] = Field(model.velocities.w; indices)
outputs[:η] = model.free_surface.η
outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

simulation.output_writers[name] = JLD2OutputWriter(model, outputs; dir,
                                                   schedule = IterationInterval(1000),
                                                   filename = output_prefix * "_fields_$rank",
                                                   with_halos = true,
                                                   overwrite_existing = true)


# MPI.Finalize()
