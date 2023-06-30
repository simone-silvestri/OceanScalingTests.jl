using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using FastSphericalHarmonics
using Dierckx

abstract type AbstractInterpolation end

struct SplineInterpolation <: AbstractInterpolation 
    passes :: Int
end

struct LinearInterpolation <: AbstractInterpolation
    passes :: Int
end

struct SpectralInterpolation{F, C} <: AbstractInterpolation
    filter_function :: F
    spectral_coeff  :: C
end

LinearInterpolation(; passes = 5) = LinearInterpolation(passes)
SpectralInterpolation(; filter_func = (l) -> exp(-l * (l+1)/ 180 / 240), spectral_coeff = nothing) =
            SpectralInterpolation(filter_func, spectral_coeff)

@kernel function _interpolate_one_level!(a, x, y, spline, ::SplineInterpolation)
    i, j = @index(Global, NTuple)
    @inbounds a[i, j] = evaluate(spline, x[i], y[j])
end

@kernel function _interpolate_one_level!(a, x, y, field, ::LinearInterpolation)
    i, j = @index(Global, NTuple)
    @inbounds a[i, j] = interpolate(field, x[i], y[j], 1)
end

function interpolate_one_level(old_array, old_grid, new_grid; interpolation_method = LinearInterpolation())
    old_field = Field{Center, Center, Nothing}(old_grid)
    set!(old_field, old_array)
    fill_halo_regions!(old_field)
    
    new_array = zeros(size(new_grid, 1), size(new_grid, 2))

    xnew = xnodes(new_grid, Center(), Center(), Center())
    ynew = ynodes(new_grid, Center(), Center(), Center())
    xold = xnodes(old_grid, Center(), Center(), Center(), with_halos = true)
    yold = ynodes(old_grid, Center(), Center(), Center(), with_halos = true)

    if interpolation_method isa SplineInterpolation
        data = parent(old_field)[:, :, 1]
        data = Spline2D(xold, yold, data)
    else
        data = old_field
    end

    launch!(CPU(), new_grid, :xy, _interpolate_one_level!, new_array, xnew, ynew, data, interpolation_method)

    return new_array
end

function interpolate_one_level_in_passes(array_old, Nxₒ, Nyₒ, Nxₙ, Nyₙ; interpolation_method = LinearInterpolation())

    # switch bathymetry to centers?

    passes = interpolation_method.passes

    ΔNx = floor((Nxₒ - Nxₙ) / passes)
    ΔNy = floor((Nyₒ - Nyₙ) / passes)

    Nx = deepcopy(Nxₒ)
    Ny = deepcopy(Nyₒ)

    @assert Nxₒ == Nxₙ + passes * ΔNx
    @assert Nyₒ == Nyₙ + passes * ΔNy
    
    array = deepcopy(array_old)
    latitude = 75.0

    for pass = 1:passes
        array_full = deepcopy(array)
        Nxₒ = Nx
        Nyₒ = Ny
        Nx -= Int(ΔNx) 
        Ny -= Int(ΔNy)
        if pass == 1
            oldlat = 89.9999999999999999
        else
            oldlat = latitude
        end
        old_grid = RectilinearGrid(size = (Nxₒ, Nyₒ), y = (-oldlat,   oldlat),   x = (-180, 180), topology = (Periodic, Bounded, Flat))
        new_grid = RectilinearGrid(size = (Nx,  Ny ), y = (-latitude, latitude), x = (-180, 180), topology = (Periodic, Bounded, Flat))
    
        @show Nxₒ, Nyₒ, Nx, Ny, pass
        array = interpolate_one_level(array_full, old_grid, new_grid; interpolation_method)
    end

    return array
end
