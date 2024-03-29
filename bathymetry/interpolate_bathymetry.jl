using PyCall
using JLD2 

"""
    interpolate_bathymetry_from_ETOPO1(resolution, maximum_latitude; 
                                       interpolation_method = LinearInterpolation(),
                                       minimum_depth = 6,
                                       etopo1_file = datadep"etopo1_bathymetry/bathymetry-ice-21600x10800.jld2")

Generate a latitude-longitude bathymetry array that spans `latitude = (-maximum_latitude, +maximum_latitude)`
with size `(360 / resolution, 2maximum_latitude / resolution)`.

Arguments
=========

- `resolution :: Float`: The lateral resolution (in degrees) of the output bathymetry.
- `latitude :: Tuple`: The south/north meridional extent of the domain.
- `longitude :: Tuple`: The west/east zonal extent of the domain.

Keyword Arguments
=================
- `interpolation_method`: Either `LinearInterpolation()` (default), SpectralInterpolation(), or SplineInterpolation().
- `minimum_depth`: The minimum depth of the bathymetry in meters. Anything less than `minimum_depth` is
                   considered land. Default `6` meters.
- `etopo1_file`: path to jld2 datafile containing the ETOPO1 bathymetry
```
"""
function interpolate_bathymetry_from_ETOPO1(resolution, latitude, longitude = (-180, 180); 
                                            interpolation_method = LinearInterpolation(),
                                            minimum_depth = 6,
                                            etopo1_file = datadep"etopo1_bathymetry/bathymetry-ice-21600x10800.jld2")

    file = jldopen(etopo1_file)
    bathy_old = Float64.(file["bathymetry"])
    
    # apparently bathymetry is reversed in the latitude direction, therefore we have to swap it
    bathy_old = reverse(bathy_old, dims = 2)

    # original file is reversed

    Nx = Int(abs(longitude[2] - longitude[1]) * resolution)
    Ny = Int(abs(latitude[2]  - latitude[1] ) * resolution)

    if interpolation_method isa SpectralInterpolation
        if interpolation_method.spectral_coeff isa nothing
            spectral_coeff = etopo1_to_spherical_harmonics(bathy_old, size(bathy_old, 2))
        else 
            spectral_coeff = interpolation_method.spectral_coeff
        end

        bathy = bathymetry_from_etopo1(Nx, Ny, spectral_coeff, interpolation_method.filter_func)
    else 
        bathy = interpolate_bathymetry_in_passes(bathy_old, size(bathy_old)..., Nx, Ny, latitude, longitude; interpolation_method)
    end

    bathy[bathy .> - minimum_depth] .= ABOVE_SEA_LEVEL

    fixed_bathymetry = remove_connected_regions(bathy)

    return fixed_bathymetry
end

function etopo1_to_spherical_harmonics(etopo1, Nmax)

    lmax = Nmax - 2

    # latitude interpolate and transpose
    etopo1_center = 0.5 *(etopo1[:, 2:end] .+ etopo1[:, 1:end-1])
    etopo1_center = etopo1_center'

    # Drop the 360 point
    etopo1_center = etopo1_center[:, 1:end-1]

    # longitude interpolate
    fft_etopo1_center = rfft(etopo1_center, 2)
    fft_etopo1_center = fft_etopo1_center[:, 1:end-1]
    etopo1_interp     = irfft(fft_etopo1_center, 2lmax +1, 2)

    # Spherical harmonic filtering
    return sph_transform(etopo1_interp)
end

function bathymetry_from_etopo1(Nφ, Nλ, spher_harm_coeff, filter)

    lmax_interp = Nφ - 1

    spher_harm_coeff_filter = zeros(lmax_interp + 1, 2lmax_interp +1)
    for l = 0:lmax_interp, m = -l:l
        spher_harm_coeff_filter[sph_mode(l, m)] = spher_harm_coeff[sph_mode(l, m)] * filter(l)
    end

    etopo1_filter = sph_evaluate(spher_harm_coeff_filter) # dimensions are Nφ, 2Nφ - 1 

    # longitude interpolate 2
    fft_etopo1_filter = rfft(etopo1_filter, 2) # dimensions are Nφ, Nφ 

    mmax_interp = Nλ ÷ 2 + 1

    fft_etopo1_interp = zeros(Complex{Float64}, Nφ, mmax_interp)

    if mmax_interp <= size(fft_etopo1_filter, 2)
        fft_etopo1_interp .= fft_etopo1_filter[:, 1:mmax_interp]
    else
        fft_etopo1_interp[:, 1:size(fft_etopo1_filter, 2)] .= fft_etopo1_filter
    end

    etopo1_final = irfft(fft_etopo1_interp, Nλ, 2)

    return etopo1_final
end

function remove_connected_regions(bat)

    bathymetry = deepcopy(bat)
    batneg     = deepcopy(bathymetry)

    batneg[batneg.>0] .= 0
    batneg[batneg.<0] .= 1

    labels = scikitimage.label(batneg, connectivity = 1)
    try
        total_elements = zeros(maximum(labels))

        for i in 1:lastindex(total_elements)
            total_elements[i] = sum(labels[labels.==i])
        end

        ocean_idx = findfirst(x -> x == maximum(x), total_elements)
        second_maximum = maximum(filter((x) -> x != total_elements[ocean_idx], total_elements))

        bering_sea_idx = findfirst(x -> x == second_maximum, total_elements)

        labels = Float64.(labels)
        labels[labels.==0] .= NaN

        for i in 1:length(total_elements)
            if (i != ocean_idx) && (i != bering_sea_idx)
                labels[labels.==i] .= NaN
            end
        end

        bathymetry .+= labels
        bathymetry[isnan.(bathymetry)] .= ABOVE_SEA_LEVEL
    catch err
        println("this is the error $err")
    end
    return bathymetry
end

write_bathymetry_to_file(resolution, bathy) = 
    JLD2.jldsave("bathymetry" * "$resolution" * ".jld2", bathymetry = bathy)
