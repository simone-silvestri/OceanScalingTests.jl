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
    @inbounds a[i, j] = interpolate(field, x[i], y[j], 0.0)
end

function interpolate_one_level(old_array, old_grid, new_grid; interpolation_method = LinearInterpolation())
    old_field = CenterField(old_grid)
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

function interpolate_bathymetry_in_passes(old_bathymetry, Nxₒ, Nyₒ, Nxₙ, Nyₙ, latitude, longitude; interpolation_method = LinearInterpolation())

    # switch bathymetry to centers?
    passes = interpolation_method.passes
    
    ΔNx = floor((Nxₒ - Nxₙ) / passes)
    ΔNy = floor((Nyₒ - Nyₙ) / passes)

    Nx_all = [Nxₒ - ΔNx * pass for pass in 1:passes-1]
    Ny_all = [Nyₒ - ΔNy * pass for pass in 1:passes-1]

    Nx_all = Int[Nxₒ, Nx_all..., Nxₙ]
    Ny_all = Int[Nyₒ, Ny_all..., Nyₙ]

    array = deepcopy(old_bathymetry)

    for pass = 2:passes+1
        array_full = deepcopy(array)
        nxₒ = Nx_all[pass-1]
        nyₒ = Ny_all[pass-1]
        nx  = Nx_all[pass]
        ny  = Ny_all[pass]
        if pass == 2
            oldlat = (-90, 90) 
            oldlon = (-180, 180)
        else
            oldlat = latitude
            oldlon = longitude
        end
	old_grid = RectilinearGrid(size = (nxₒ, nyₒ, 1), y = (  oldlat[1],   oldlat[2]), x = (   oldlon[1],    oldlon[2]), z = (-1, 1), topology = (Periodic, Bounded, Bounded))
	new_grid = RectilinearGrid(size = (nx,  ny , 1), y = (latitude[1], latitude[2]), x = (longitude[1], longitude[2]), z = (-1, 1), topology = (Periodic, Bounded, Bounded))
    
        @show nxₒ, nyₒ, nx, ny, pass
        array = interpolate_one_level(array_full, old_grid, new_grid; interpolation_method)
    end

    return array
end
