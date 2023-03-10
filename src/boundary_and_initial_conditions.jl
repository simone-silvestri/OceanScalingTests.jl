using Oceananigans.Units

initialize_model!(model, ::Val{:Quiescent})   = nothing

function initialize_model!(model, ::Val{:DoubleDrake})
    Depth = model.grid.Lz
    exp_length = Depth / 15.0

    @inline init_temperature(λ, φ, z) = exponential_profile(z; Lz = Depth, h = exp_length) * max(T_reference(φ), 2.0)

    set!(model, T = init_temperature, S = 35.0)
end

function initialize_model!(model, ::Val{:RealisticOcean})
    grid = model.grid

    rx = grid.architecture.local_rank
    nx = size(grid, 1) 

    T_init = jldopen("../data/T_init.jld2")[1+rx*nx:(rx+1)*nx, :]
    S_init = jldopen("../data/S_init.jld2")[1+rx*nx:(rx+1)*nx, :]
    
    set!(model, T = T_init, S = S_init)
end

function set_boundary_conditions(::Val{:Quiescent}, Ly)

    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function set_boundary_conditions(::Val{:DoubleDrake}, Ly)

    S_reference = 35.0 # reference salinity for salinity flux (psu)
    u_coeffs = wind_stress_coefficients(Ly/2)
    S_coeffs = salinity_flux_coefficients(Ly/2, S_reference)

    λ = 5meters / 7days # Pumping velocity of temperature restoring (ms⁻¹)
    u_top_bc = FluxBoundaryCondition(surface_stress_x,      discrete_form=true, parameters=u_coeffs)
    S_top_bc = FluxBoundaryCondition(surface_salinity_flux, discrete_form=true, parameters=S_coeffs)
    T_top_bc = FluxBoundaryCondition(T_relaxation,          discrete_form=true, parameters=λ)

    μ = 0.003 # linear drag coefficient (ms⁻¹)
    u_bot_bc = FluxBoundaryCondition(u_linear_bottom_drag, discrete_form=true, parameters=μ)
    v_bot_bc = FluxBoundaryCondition(v_linear_bottom_drag, discrete_form=true, parameters=μ)

    T_bcs = FieldBoundaryConditions(top=T_top_bc)
    S_bcs = FieldBoundaryConditions(top=S_top_bc)
    u_bcs = FieldBoundaryConditions(bottom=u_bot_bc, top=u_top_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bot_bc)

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function set_boundary_conditions(::Val{:RealisticOcean}, Ly)
    # This function assumes data is available to load from the folder `../data`
    τˣ = nothing
    τʸ = nothing
    Tˢ = nothing
    Sˢ = nothing

    μ = 0.001 # Quadratic drag coefficient (ms⁻¹)
    u_bot_bc = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=μ)
    v_bot_bc = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_bot_bc = FluxBoundaryCondition(u_immersed_quadratic_bottom_drag, discrete_form=true, parameters=μ)
    v_immersed_bot_bc = FluxBoundaryCondition(v_immersed_quadratic_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = u_immersed_bot_bc)
    u_immersed_bc = ImmersedBoundaryCondition(bottom = v_immersed_bot_bc)

    T_bcs = FieldBoundaryConditions(top=T_top_bc)
    S_bcs = FieldBoundaryConditions(top=S_top_bc)
    u_bcs = FieldBoundaryConditions(bottom=u_bot_bc, top=u_top_bc, immersed = u_immersed_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bot_bc, top=v_top_bc, immersed = v_immersed_bc)

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end