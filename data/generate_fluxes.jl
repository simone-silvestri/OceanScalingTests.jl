using Statistics: dot
using DataDeps
using NetCDF
using JLD2 

include("interpolation_utils.jl")

iters = Int[]
final_year  = parse(Int, get(ENV, "FINALYEAR", "2000"))
final_month = parse(Int, get(ENV, "FINALMONTH", "12")) 
for year in 1995:final_year-1
    for month in 1:12
        for day in OceanScalingTests.monthly_days(year)[month]
            push!(iters, parse(Int, string(year) * string(month, pad=2) * string(day, pad=2)))
        end
    end
end
year = final_year
for month in 1:final_month
    for day in OceanScalingTests.monthly_days(year)[month]
        push!(iters, parse(Int, string(year) * string(month, pad=2) * string(day, pad=2)))
    end
end

# Starting from the 1st of January 1992
it_collection = []
for i = 1:length(iters)÷5
    push!(it_collection, iters[1+(i-1)*5:i*5+1])
end

jldsave("list_of_days.jld2", list = it_collection)

YField = Field{<:Center, <:Face, <:Center}
set_field!(f::Field, a) = set!(f, a)

function set_field!(f::YField, a)
    d = zeros(size(f)...)
    d[:, 1:end-1] .= a

    set!(f, d)

    return nothing
end

function read_and_interpolate_quarter_flux(name, iterations, Nx, Ny, Nj = 1; arch = GPU(), 
                                           max_val = nothing, location = (Center, Center, Center))

    grid_cs510 = RectilinearGrid(arch, size = (1440, 720, 1), y = (-90, 90), x = (-180, 180), z = (0, 1), topology = (Periodic, Bounded, Bounded)) 
    cs_field = Field{location...}(grid_cs510)

    interp = if location[2] == Face
        zeros(Nx, Ny+1, 6)
    else
        zeros(Nx, Ny, 6)
    end

    for j in 1:Nj

        Δφ = 150 / Nj

        latitude = (- 75 + Δφ * (j - 1), - 75 + Δφ * j)
        
        Nφ = Int(Ny / Nj)
                    
        new_grid = RectilinearGrid(arch; size = (Nx, Nφ, 1), y = latitude, x = (-180, 180), z = (0, 1), topology = (Periodic, Bounded, Bounded)) 
        my_field = Field{location...}(new_grid)
        my_tmp2  = Field{location...}(new_grid)

        tmp_interp = zeros(size(my_field)[1:2]..., length(iterations))

        for (idx, iter) in enumerate(iterations)

            file = "$(name).1440x720.$(iter).nc"

            # Downloading files
            url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/$(name)_daily.nc/$(file)"
            cmd = `wget --http-user=ssilvestri --http-passwd=ZZjQeLy7oIHwvqMWvM8y $(url)`
                            
            !isfile(file) && run(cmd)

            @info "interpolating file $file"

            set_field!(cs_field, ncread(file, name))
            fix_max_val!(cs_field, max_val)

            fill_halo_regions!(cs_field)

            for step in 1:50
                propagate_field!(cs_field)
            end

            horizontal_interpolate!(my_field, cs_field)

            fill_halo_regions!(my_field)
            
            for step in 1:5
                horizontal_filter!(my_tmp2, my_field)
                horizontal_filter!(my_field, my_tmp2)
            end

            tmp_interp[:, :, idx] .= Array(interior(my_field))
        end

        jrange = UnitRange(1 + (j - 1) * Nφ, size(my_field, 2) + (j - 1) * Nφ)
        interp[:, jrange, :] .= tmp_interp
    end
    return interp
end

function transpose_flux!(var, tmp)
    nx = size(tmp, 1)
    tmp .= var[nx+1:end, :, :]
    var[nx+1:end, :, :] .= var[1:nx, :, :]
    var[1:nx, :, :]   .= tmp
end

function generate_fluxes(resolution; arch = GPU())

    # grid size
    Nx = Int(360 * resolution) 
    Ny = Int(150 * resolution)

    tmp  = zeros(Nx÷2, Ny,   6)
    tmpy = zeros(Nx÷2, Ny+1, 6)

    for (idx, iterations) in enumerate(it_collection)

        τx = read_and_interpolate_quarter_flux("oceTAUX",  iterations, Nx, Ny; arch, max_val = 1e8, location = (Face, Center, Center))
        τy = read_and_interpolate_quarter_flux("oceTAUY",  iterations, Nx, Ny; arch, max_val = 1e8, location = (Center, Face, Center))
        Fs = read_and_interpolate_quarter_flux("oceFWflx", iterations, Nx, Ny; arch, max_val = 1e8)
        Qs = read_and_interpolate_quarter_flux("oceQnet",  iterations, Nx, Ny; arch, max_val = 1e8)

        @info "Correcting fluxes"
        Qs .*= (1 / 1000 / 3991)
        Fs .*= (1 / 1000 * 35)

        Qs .= -Qs
        τx .= -τx ./ 1000
        τy .= -τy ./ 1000

        transpose_flux!(τx, tmp)
        transpose_flux!(Qs, tmp)
        transpose_flux!(Fs, tmp)
        transpose_flux!(τy, tmpy)

        @info "saving down fluxes $idx"
        jldsave("fluxes_$(idx).jld2", τx = τx, τy = τy, Qs = Qs, Fs = Fs)
    end
end
