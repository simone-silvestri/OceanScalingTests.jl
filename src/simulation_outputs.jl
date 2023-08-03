using Oceananigans.Units
using Oceananigans.Grids: halo_size
using Oceananigans.BoundaryConditions
using Oceananigans.OutputWriters: set_time_stepper!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField

import Oceananigans.OutputWriters: set!, set_time_stepper_tendencies!

function set_outputs!(simulation, ::Val{Experiment}; overwrite_existing = true, surface_time = 1days, checkpoint_time = 5days, output_dir = "./") where Experiment

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
                                                           dir = output_dir, 
                                                           schedule = TimeInterval(surface_time),
                                                           filename = "$(Experiment)_fields_$rank",
                                                           with_halos = true,
                                                           overwrite_existing)

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            dir = output_dir, 
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
                @info name size(parent_data) size(model_field)
		copyto!(model_field.data.parent, parent_data)
            end
        end

        load_free_surface!(model_fields.η, size(model_fields.T), file)
        set_time_stepper!(model.timestepper, file, model_fields)

        checkpointed_clock = file["clock"]

        # Update model clock
        model.clock.iteration = checkpointed_clock.iteration
        model.clock.time = checkpointed_clock.time
    end

    return nothing
end

function load_free_surface!(η, N, file)
    Sη = size(file["η/data"])[1:2]
    N  = N[1:2] 

    Hx, Hy = Int.((Sη .- N) ./ 2)
    @show Hx Hy 
    η_file = file["η/data"][Hx+1:end-Hx, Hy+1:end-Hy, :]
    set!(η, η_file)
    fill_halo_regions!(η)

    return nothing
end

function set_time_stepper_tendencies!(timestepper, file, model_fields)
    for name in (:u, :v, :T, :S)
        if string(name) ∈ keys(file["timestepper/Gⁿ"]) # Test if variable tendencies exist in checkpoint
            # Tendency "n"
            parent_data = file["timestepper/Gⁿ/$name/data"]

            tendencyⁿ_field = timestepper.Gⁿ[name]
            copyto!(tendencyⁿ_field.data.parent, parent_data)

            # Tendency "n-1"
            parent_data = file["timestepper/G⁻/$name/data"]

            tendency⁻_field = timestepper.G⁻[name]
            copyto!(tendency⁻_field.data.parent, parent_data)
        else
            @warn "Tendencies for $name do not exist in checkpoint and could not be restored."
        end
    end

    return nothing
end
