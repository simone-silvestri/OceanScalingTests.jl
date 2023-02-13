using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

function set_outputs!(simulation, ::Val{:DoubleDrake})

    model   = simulation.model

    outputs = Dict()
    indices = (:, :, model.grid.Nz) 

    outputs[:u] = Field(model.velocities.u; indices)
    outputs[:v] = Field(model.velocities.v; indices)
    outputs[:w] = Field(model.velocities.w; indices)
    outputs[:T] = Field(model.tracers.T;    indices)
    outputs[:S] = Field(model.tracers.S;    indices)
    outputs[:η] = model.free_surface.η
    outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

    rank = model.grid.architecture.local_rank

    simulation.output_writers[:surface] = JLD2OutputWriter(model, outputs;
                                                           schedule = IterationInterval(1000),
                                                           filename = "double_drake_fields_$rank",
                                                           with_halos = true,
                                                           overwrite_existing = true)

    simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                            schedule = IterationInterval(500000),
                                                            prefix = "double_drake_checkpoint_$rank",
                                                            overwrite_existing = true)
end

set_outputs!(simulation, ::Val{:Quiescent}) = nothing