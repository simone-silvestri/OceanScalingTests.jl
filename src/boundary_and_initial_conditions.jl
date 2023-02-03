
initialize_model!(model, ::Val{:Quiescent})   = nothing

function initialize_model!(model, ::Val{:DoubleDrake})
    Depth = model.grid.Lz
    exp_length = Depth / 5.0

    @inline init_temperature(λ, φ, z) = initial_temperature(λ, φ, z; Lz = Depth, h = exp_length)

    set!(model, T = init_temperature, S = 35.0)
end

function boundary_conditions(::Val{:Quiescent})

    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function boundary_conditions(::Val{:DoubleDrake})

    u_coeffs = wind_stress_coefficients(80.0)
    S_coeffs = salinity_flux_coefficients(80.0, 35.0)

    u_top_bc = FluxBoundaryCondition(surface_stress_x,      discrete_form=true, parameters=(; coeffs = u_coeffs))
    S_top_bc = FluxBoundaryCondition(surface_salinity_flux, discrete_form=true, parameters=(; coeffs = S_coeffs))
    T_top_bc = FluxBoundaryCondition(T_relaxation,          discrete_form=true, parameters=30days)

    u_bot_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form=true, parameters = (; μ = 0.1))
    v_bot_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form=true, parameters = (; μ = 0.1))

    T_bcs = FieldBoundaryConditions(top=T_top_bc)
    S_bcs = FieldBoundaryConditions(top=S_top_bc)
    u_bcs = FieldBoundaryConditions(bottom=u_bot_bc, top=u_top_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bot_bc)

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end
