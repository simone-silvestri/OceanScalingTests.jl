using Oceananigans.Grids: ynode

function linear_z_faces(Nz, Depth; first_Δ = 5)

    m = (Depth - (Nz * first_Δ)) / sum(1:Nz-1)
    Δ = m .* (0:Nz-1) .+ first_Δ

    z_faces = zeros(Nz+1)
    for k in Nz:-1:1
        z_faces[k] = z_faces[k+1] - Δ[Nz - k + 1]
    end

    return z_faces
end

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_z_faces(Nz, Depth; h = Nz / 3)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - Depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
end

function double_drake_bathymetry(λ, φ) 
    if φ > -35
        (λ >  0 && λ < 1)  && return 0.0
        (λ > 90 && λ < 91) && return 0.0
    end

    return -10000.0
end

# Assumes there is a jld2 file called bathymetry.jld2 in the data folder
function realistic_bathymetry(grid)
    rx = grid.architecture.local_rank
    nx = size(grid, 1) 
    return jldopen("../data/bathymetry.jld2")[1+rx*nx:(rx+1)*nx, :]
end

@inline function cubic_profile(x1, x2, y1, y2, d1, d2)
    A = [ x1^3 x1^2 x1 1.0
          x2^3 x2^2 x2 1.0
          3*x1^2 2*x1 1.0 0.0
          3*x2^2 2*x2 1.0 0.0]
          
    b = [y1, y2, d1, d2]

    coeff = A \ b

    return tuple(coeff...)
end

# Coefficients for a piecewise cubic interpolation of wind stress 
@inline function wind_stress_coefficients(south_north_limit)
    below_45_coeffs = cubic_profile(-south_north_limit, -45.0, 0.0, 0.2, 0.0, 0.0) ./ 1000
    below_15_coeffs = cubic_profile(-45.0, -15.0, 0.2, -0.1, 0.0, 0.0) ./ 1000
    below_00_coeffs = cubic_profile(-15.0, 0.0, -0.1, -0.02, 0.0, 0.0) ./ 1000
    above_00_coeffs = cubic_profile(0.0, 15.0, -0.02, -0.1, 0.0, 0.0) ./ 1000
    above_15_coeffs = cubic_profile(15.0, 45.0, -0.1, 0.1, 0.0, 0.0) ./ 1000
    above_45_coeffs = cubic_profile(45.0, south_north_limit, 0.1, 0.0, 0.0, 0.0) ./ 1000
    
    return (below_45_coeffs, below_15_coeffs, below_00_coeffs, above_00_coeffs, above_15_coeffs, above_45_coeffs)
end

# Coefficients for a piecewise cubic interpolation of surface salinity flux
@inline function salinity_flux_coefficients(south_north_limit, initial_salinity)
    below_20_coeffs = cubic_profile(-south_north_limit, -20.0, -2e-8, 2e-8, 0.0, 0.0) .* initial_salinity
    below_00_coeffs = cubic_profile(-20.0, 0.0, 2e-8, -4e-8, 0.0, 0.0) .* initial_salinity
    above_00_coeffs = cubic_profile(0.0, 20.0, -4e-8, 2e-8, 0.0, 0.0) .* initial_salinity
    above_20_coeffs = cubic_profile(20.0, south_north_limit, 2e-8, -2e-8, 0.0, 0.0) .* initial_salinity

    return (below_20_coeffs, below_00_coeffs, above_00_coeffs, above_20_coeffs)
end
    
@inline function wind_stress(φ, coeffs) 
    coeff = ifelse(φ < -45, coeffs[1], 
            ifelse(φ < -15, coeffs[2], 
            ifelse(φ <   0, coeffs[3], 
            ifelse(φ <  15, coeffs[4], 
            ifelse(φ <  45, coeffs[5], 
                            coeffs[6]))))) 
    
    return coeff[1] * φ^3 + coeff[2] * φ^2 + coeff[3] * φ + coeff[4]
end

@inline function salinity_flux(φ, coeffs) 
    coeff = ifelse(φ < -20, coeffs[1], 
            ifelse(φ <   0, coeffs[2], 
            ifelse(φ <  20, coeffs[3], 
                            coeffs[4])))
    
    return coeff[1] * φ^3 + coeff[2] * φ^2 + coeff[3] * φ + coeff[4]
end

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

@inline u_linear_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1]
@inline v_linear_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1]

@inline function surface_stress_x(i, j, grid, clock, fields, p)
    φ = ynode(Center(), j, grid)
    return wind_stress(φ, p)
end

@inline function surface_salinity_flux(i, j, grid, clock, fields, p)
    φ = ynode(Center(), j, grid)
    return salinity_flux(φ, p)
end

@inline T_reference(φ) = max(0.0, 30.0 * cos(1.2 * π * φ / 180))

@inline function T_relaxation(i, j, grid, clock, fields, λ)
    φ = ynode(Center(), j, grid)
    return @inbounds λ * (fields.T[i, j, grid.Nz] - T_reference(φ))
end

# Fluxes are saved as [Nx, Ny, Nt] where Nt = 1:10 and represents day 0 to day 9
@inline function flux_interpolate_array(i, j, grid, clock, fields, p)
    n  = mod(clock.time, 9days) + 1
    if n == 1
        return p[i, j, n]
    end
    n₁ = Int(floor(n))
    n₂ = n₁ + 1
    return p[i, j, n₁] * (n₂ - n) + p[i, j, n₂] * (n - n₁)
end

# Update fluxes through a Callback every 9days
@inline function update_fluxes(sim)

    filenum = Int(sim.model.clock.time ÷ 9days)
    file    = jldopen("data/fluxes_$(filenum).jld2")

    model = sim.model
    grid  = model.grid

    u_top = model.velocities.u.boundary_conditions.top.condition 
    v_top = model.velocities.u.boundary_conditions.top.condition 
    T_top = model.velocities.u.boundary_conditions.top.condition 
    S_top = model.velocities.u.boundary_conditions.top.condition 

    rx = grid.architecture.local_rank
    nx = size(grid, 1) 

    u_top.parameters .= file["τx"][1+rx*nx:(rx+1)*nx, :, :]
    v_top.parameters .= file["τy"][1+rx*nx:(rx+1)*nx, :, :]

    T_top.parameters .= file["Q"][1+rx*nx:(rx+1)*nx, :, :]
    S_top.parameters .= file["F"][1+rx*nx:(rx+1)*nx, :, :]

    return nothing
end

