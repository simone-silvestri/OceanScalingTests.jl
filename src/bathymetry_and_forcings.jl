using Oceananigans.Grids: φnode
using Oceananigans.Distributed
using Oceananigans.Distributed: DistributedGrid, partition_global_array
using Oceananigans.Operators
using DataDeps

function z_from_ecco(Nz, Depth)

    if Nz != 48
	    throw(ArgumentError("Not ecco grid!!"))
    end
    path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

    dh = DataDep("quarter_degree_near_global_lat_lon",
      "Forcing data for global latitude longitude simulation",
       path * "z_faces-50-levels.jld2"
    )

    DataDeps.register(dh)

    datadep"quarter_degree_near_global_lat_lon"
    
    datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/z_faces-50-levels.jld2"
    file_z_faces = jldopen(datadep_path)
    
    # Stretched faces taken from ECCO Version 4 (50 levels in the vertical)
    return file_z_faces["z_faces"][3:end];
end

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

function exponential_z_faces(Nz, Depth; h = Nz / 4.5)

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
realistic_bathymetry(grid::DistributedGrid, resolution) = 
      partition_global_array(architecture(grid), eltype(grid).(jldopen("bathymetry/bathymetry$(resolution).jld2")["bathymetry"]), size(grid))

realistic_bathymetry(grid, resolution) = eltype(grid).(jldopen("bathymetry/bathymetry$(resolution).jld2")["bathymetry"])

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
@inline spᶜᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(ℑxᶜᵃᵃ(i, j, k, grid, ϕ², Φ.u) + ℑyᵃᶜᵃ(i, j, k, grid, ϕ², Φ.v))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

@inline u_linear_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1]
@inline v_linear_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1]

@inline function surface_stress_x(i, j, grid, clock, fields, p)
    φ = φnode(j, grid, Center())
    return wind_stress(φ, p)
end

@inline function surface_salinity_flux(i, j, grid, clock, fields, p)
    φ = φnode(j, grid, Center())
    return salinity_flux(φ, p)
end

@inline T_reference(φ) = max(0.0, 30.0 * cos(1.2 * π * φ / 180))

@inline function T_relaxation(i, j, grid, clock, fields, λ)
    φ = φnode(j, grid, Center())
    return @inbounds λ * (fields.T[i, j, grid.Nz] - T_reference(φ))
end

# Fluxes are saved as [Nt, Nx, Ny] where Nt = 1:6 and represents day 0 to day 5
@inline function flux_from_interpolated_array(i, j, grid, clock, fields, p)
    time_in_days = clock.time / 1days
    n  = mod(time_in_days, 5) + 1
    n₁ = Int(floor(n))
    n₂ = Int(n₁ + 1)    
    return p[i, j, n₁] * (n₂ - n) + p[i, j, n₂] * (n - n₁)
end
