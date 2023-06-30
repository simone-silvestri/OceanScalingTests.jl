using Oceananigans 
using Oceananigans.BoundaryConditions
using Oceananigans.Utils: launch!
using Oceananigans.Architectures
using Oceananigans.Architectures: architecture
using Oceananigans.Distributed
using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll
using Statistics: dot
using DataDeps
using NetCDF
using JLD2 

@kernel function _propagate_field!(field, Nx, Ny)
    i, j, k = @index(Global, NTuple)
    nw = ifelse(i == 1 , field[Nx, j, k], field[i - 1, j, k])
    ns = ifelse(j == 1 , 0.0,             field[i, j - 1, k])
    ne = ifelse(i == Nx, field[1, j, k],  field[i + 1, j, k])
    nn = ifelse(j == Ny, 0.0,             field[i, j + 1, k])
    nn = (nw, ne, nn, ns)
    pos = Int.(nn .!= 0)

    if (field[i, j, k] == 0) & (sum(pos) > 0)
        field[i, j, k] = dot(pos, nn) / sum(pos) 
    end
end

propagate_field!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _propagate_field!, field, Nx, Ny)

@kernel function _extend_vertically!(field, Nz)
    i, j = @index(Global, NTuple)

    @unroll for k in Nz-1:-1:1
        if field[i, j, k] == 0.0
            field[i, j, k] = field[i, j, k+1]
        end
    end
end

extend_vertically!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _extend_vertically!, field, size(grid, 3))

@kernel function _horizontal_filter!(new_field, field)
    i, j, k = @index(Global, NTuple)

    new_field[i, j, k] = field[i, j, k]
    nn = (field[i, j, k], field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k])
    non_null  = Int.(nn .!= 0)

    if sum(non_null) > 0
        new_field[i, j, k] = sum(nn) / sum(non_null)
    end
end

horizontal_filter!(new_field, field) =
    launch!(architecture(field.grid), field.grid, :xyz, _horizontal_filter!, new_field, field)

@kernel function _fix_max_val!(field, max_val)
    i, j, k = @index(Global, NTuple)

    if abs(field[i, j, k]) > max_val
        field[i, j, k] = 0.0
    end
end

fix_max_val!(field, ::Nothing) = nothing

fix_max_val!(field, max_val) =
    launch!(architecture(field.grid), field.grid, :xyz, _fix_max_val!, field, max_val)

using Oceananigans.Fields: interpolate, location, instantiated_location
using Oceananigans.Grids: xnode, ynode, znode

@kernel function _horizontal_interpolate!(new_field, old_field, new_grid, old_grid, loc)
    i, j, k = @index(Global, NTuple)
    x = xnode(loc[1], i, new_grid)
    y = ynode(loc[2], j, new_grid)
    z = znode(loc[3], k, old_grid)
    
    new_field[i, j, k] = interpolate(old_field, loc..., old_grid, x, y, z)
end

@inline function horizontal_interpolate!(new_field, old_field)
    
    new_grid = new_field.grid  
    old_grid = old_field.grid  

    arch     = architecture(new_grid)
    location = instantiated_location(new_field)

    launch!(arch, new_grid, :xyz, _horizontal_interpolate!, new_field, old_field, new_grid, old_grid, location)
    
    return nothing
end

# Doing the year 1995-2020
leap_year_days(year) = year == 1996 || 
                       year == 2000 || 
                       year == 2004 || 
                       year == 2008 || 
                       year == 2012 || 
                       year == 2016 || 
                       year == 2020 ? UnitRange(1, 29) : UnitRange(1, 28)

monthly_days(year) = [1:31, leap_year_days(year), 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

iters = Int[]
final_year = parse(Int, get(ENV, "FINALYEAR", "2000"))

for year in 1995:final_year
    for month in 1:12
        for day in monthly_days(year)[month]
            push!(iters, parse(Int, string(year) * string(month, pad=2) * string(day, pad=2)))
        end
    end
end

# Starting from the 1st of January 1992
it_collection = []
for i = 1:length(iters)÷10
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

    grid_cs510 = LatitudeLongitudeGrid(arch, size = (1440, 720, 1), latitude = (-90, 90), longitude = (-180, 180), z = (0, 1)) 
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
                    
        new_grid = LatitudeLongitudeGrid(arch; size = (Nx, Nφ, 1), latitude, longitude = (-180, 180), z = (0, 1)) 
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
