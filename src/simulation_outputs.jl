using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.Units
using Oceananigans.OutputWriters: set_time_stepper!
using Oceananigans.BoundaryConditions

import Oceananigans.OutputWriters: set!

function set_outputs!(simulation, ::Val{Experiment}; overwrite_existing = true, surface_time = 1days, checkpoint_time = 5days) where Experiment

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
                                                           schedule = TimeInterval(surface_time),
                                                           filename = "$(Experiment)_fields_$rank",
                                                           with_halos = true,
                                                           overwrite_existing)

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = "$(Experiment)_checkpoint_$rank",
                                                            overwrite_existing)
end

set_outputs!(simulation, ::Val{:Quiescent}) = nothing

function set!(model, filepath::AbstractString)

    jldopen(filepath, "r") do file

        # Do NOT validate the grid!
        model_fields = prognostic_fields(model)

        for name in (:u, :v, :T, :S)
            if string(name) ∈ keys(file) # Test if variable exist in checkpoint.
                model_field = model_fields[name]
                parent_data = file["$name/data"] #  Allow different halo size by loading only the interior
                copyto!(model_field.data.parent, parent_data)
            end
        end

        set_time_stepper!(model.timestepper, file, model_fields)

        η      = model_fields.η
        η_file = file["η/data"][74:end-73, 8:end-7, :]
        set!(η, η_file)
        fill_halo_regions!(η)

        checkpointed_clock = file["clock"]

        # Update model clock
        model.clock.iteration = checkpointed_clock.iteration
        model.clock.time = checkpointed_clock.time
    end

    return nothing
end
