using Statistics: dot
using DataDeps
using NetCDF
using JLD2 

include("interpolation_utils.jl")

iters = Int[]
final_year  = parse(Int, get(ENV, "FINALYEAR", "1996"))
final_month = parse(Int, get(ENV, "FINALMONTH", "1")) 
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

iters_restoring = [iters[i] for i in 1:3:length(iters)]

# Starting from the 1st of January 1992
it_collection = []
for i = 1:length(iters)÷5
    push!(it_collection, iters[1+(i-1)*5:i*5+1])
end

# Starting from the 1st of January 1992
it_collection_restoring = []
for i = 1:length(iters_restoring)÷5
    push!(it_collection_restoring, iters_restoring[1+(i-1)*5:i*5+1])
end

jldsave("list_of_days.jld2", list_fluxes = it_collection, list_restoring = it_collection_restoring)

YField = Field{<:Center, <:Face, <:Center}
set_field!(f::Field, a) = set!(f, a)

function set_field!(f::YField, a)
    d = zeros(size(f)...)
    d[:, 1:end-1] .= a

    set!(f, d)

    return nothing
end

function read_and_interpolate_quarter_flux(name, iterations, Nx, Ny, Nj = 1; arch = GPU(), 
                                           max_val = nothing, location = (Center, Center, Center),
                                           full_field = false)

    grid_cs510 = RectilinearGrid(arch, size = (1440, 720, 1), y = (-90, 90), x = (-180, 180), z = (0, 1), topology = (Periodic, Bounded, Bounded)) 
    cs_field = Field{location...}(grid_cs510)
    cs_tmp   = Field{location...}(grid_cs510)

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

            file   = full_field ? "$(name).1440x720x50.$(iter).nc" : "$(name).1440x720.$(iter).nc"
            suffix = full_field ? "" : "_daily"
            # Downloading files
            url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/$(name)$(suffix).nc/$(file)"
            cmd = `wget --http-user=ssilvestri --http-passwd=ZZjQeLy7oIHwvqMWvM8y $(url)`
                            
            !isfile(file) && run(cmd)

            @info "interpolating file $file"

            set_field!(cs_field, ncread(file, name)[:, :, 1])
            fix_max_val!(cs_field, max_val)
            fill_halo_regions!(cs_field)

            propagate_field!(cs_field, cs_tmp)
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
	    if !isfile("fluxes_$(idx).jld2")
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
        else
            @info "fluxes_$(idx).jld2 already exists!!"
        end
    end
end

function generate_restoring(resolution; arch = GPU())

    # grid size
    Nx = Int(360 * resolution) 
    Ny = Int(150 * resolution)

    tmp  = zeros(Nx÷2, Ny, 6)

    for (idx, iterations) in enumerate(it_collection_restoring)
	    if !isfile("restoring_$(idx).jld2")
            Tr = read_and_interpolate_quarter_flux("THETA", iterations, Nx, Ny; full_field = true, arch, max_val = 1e8)
            Sr = read_and_interpolate_quarter_flux("SALT",  iterations, Nx, Ny; full_field = true, arch, max_val = 1e8)

            @info "Correcting restoring"
            transpose_flux!(Tr, tmp)
            transpose_flux!(Sr, tmp)

            @info "saving down restoring $idx"
            jldsave("restoring_$(idx).jld2", Tr = Tr, Sr = Sr)
        else
            @info "restoring_$(idx).jld2 already exists!!"
        end
    end
end
