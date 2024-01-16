using FFTW

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

import Base

Base.:(+)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .+ t.spec, s.freq)
Base.:(*)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .* t.spec, s.freq)
Base.:(/)(s::Spectrum, t::Int)      = Spectrum(s.spec ./ t, s.freq)

Base.real(s::Spectrum) = Spectrum(real.(s.spec), s.freq)
Base.abs(s::Spectrum) = Spectrum(abs.(s.spec), s.freq)

@inline onefunc(args...)  = 1.0
@inline hann_window(n, N) = sin(π * n / N)^2 

function power_spectrum_1d_x(var, dx; windowing = onefunc)

    Nx = length(x)
    Nfx = Int64(Nx)
    
    spectra = zeros(ComplexF64, Int(Nfx/2))
    
    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    windowed_var = [var[i] * windowing(i, Nfx) for i in 1:Nfx]
    fourier      = fft(windowed_var) / Nfx
    spectra[1]  += fourier[1] .* conj(fourier[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += 2.0 * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
    end
    return Spectrum(spectra, freqs)
end

function _compute_zonal_spectra!(Uspec, Vspec, ηspec, dx, irange, jrange, krange, windowing, u, v, η)
    @inbounds for (j′, j) in enumerate(jrange)
        @inbounds for (k′, k) in enumerate(krange)
            Uspec[j′, k′]  = power_spectrum_1d_x(u[irange, j, k], dx; windowing)
            Vspec[j′, k′]  = power_spectrum_1d_x(v[irange, j, k], dx; windowing)
        end
        ηspec[j′]  = power_spectrum_1d_x(η[:, j, 101], x; windowing)
    end
end

function _update_spectra!(Ufinal, Vfinal, ηfinal, Uspec, Vspec, ηspec)
    @inbounds for j in 1:length(Ufinal[:, 1])
        @inbounds for k in 1:length(Ufinal[1, :])
            Ufinal[j, k].spec  .+= Uspec[j, k].spec
            Vfinal[j, k].spec  .+= Vspec[j, k].spec
        end
        ηfinal[j].spec  .+= ηspec[j].spec  
    end
end

function compute_spectra(filenames, irange, jrange, krange; Δx = 1/12, windowing = hann_window)

    Ny = length(jrange)
    Nz = length(krange)

    U  = Array{Spectrum}(undef, Ny, Nz)
    V  = Array{Spectrum}(undef, Ny, Nz)
    Ω  = Array{Spectrum}(undef, Ny, Nz)
    η  = Array{Spectrum}(undef, Ny)
    
    Uspec  = Array{Spectrum}(undef, Ny, Nz)
    Vspec  = Array{Spectrum}(undef, Ny, Nz)
    ηspec  = Array{Spectrum}(undef, Ny)
    
    file = jldopen(filenames[1])

    u = file["u"]
    v = file["v"]
    η = file["η"]

    _compute_zonal_spectra!(U, V, η, Δx, irange, jrange, krange, windowing, u, v, η)

    for filename in filenames[2:end]
        @info "doing file $filename"

        file = jldopen(filenames[1])

        u = file["u"]
        v = file["v"]
        η = file["η"]

        _compute_zonal_spectra!(Uspec, Vspec, ηspec, Δx, irange, jrange, krange, windowing, u, v, η)
        _update_spectra!(U, V, η, Uspec, Vspec, ηspec)
    end

    @info "finished"
    
    return (; U, V, B, Ω, WB, ST, PV)
end

using JLD2

myfiles = readdir("./")
myfiles = filter(x -> length(x) >= 21, myfiles)
myfiles = filter(x -> x[1:21] == "compressed_iteration_", myfiles)

@inline remove_last_character(s) = s[1:end-1]
myfiles = remove_last_character.(myfiles)

numbers = parse.(Int, filter.(isdigit, myfiles))
numbers = sort(numbers)

monthly_days(year) = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

total_days = [sum([monthly_days(1995)[i][end] for i in 1:j]) for j in 1:12]

months    = zeros(12)
filenames = Tuple("compressed_iteration_" * string(num) * ".jld2" for num in numbers)

Δx = 1/12

x_range_southern_ocean = 1:4320 
x_range_kuroshio       = 3870:4320
x_range_gulf_stream    = 1250:1900

j_range_southern_ocean = 150:210
j_range_kuroshio       = 1200:1400
j_range_gulf_stream    = 1200:1370

k_range_kuroshio       = 90:100
k_range_southern_ocean = 90:100
k_range_gulf_stream    = 90:100

southern_ocean_spectra = compute_spectra(filenames[15:end], x_range_southern_ocean, 
                                                            j_range_southern_ocean, 
                                                            k_range_southern_ocean)
jldsave("kuroshio_spectra.jld2", spectra = southern_ocean_spectra)

kuroshio_spectra = compute_spectra(filenames[15:end], x_range_kuroshio, 
                                                      j_range_kuroshio, 
                                                      k_range_kuroshio; 
                                                      windowing = hann_window)
jldsave("kuroshio_spectra.jld2", spectra = kuroshio_spectra)

gulf_stream_spectra = compute_spectra(filenames[15:end], x_range_gulf_stream, 
                                                         j_range_gulf_stream, 
                                                         k_range_gulf_stream; 
                                                         windowing = hann_window)
jldsave("kuroshio_spectra.jld2", spectra = gulf_stream_spectra)