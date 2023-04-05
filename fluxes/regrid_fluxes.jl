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

arch = GPU()

resolution = parse(Int, get(ENV, "RESOLUTION", "3"))

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

function propagate_field!(field) 
    grid = field.grid
    arch = architecture(grid)

    Nx, Ny, _ = size(grid)

    event = launch!(arch, grid, :xyz, _propagate_field!, field, Nx, Ny)
    wait(Oceananigans.Architectures.device(arch), event)

    return nothing
end

@kernel function _extend_vertically!(field, Nz)
    i, j = @index(Global, NTuple)

    @unroll for k in Nz-1:-1:1
        if field[i, j, k] == 0.0
            field[i, j, k] = field[i, j, k+1]
        end
    end
end

function extend_vertically!(field) 
    grid = field.grid
    arch = architecture(grid)

    event = launch!(arch, grid, :xy, _extend_vertically!, field, size(grid, 3))
    wait(Oceananigans.Architectures.device(arch), event)

    return nothing
end

@kernel function _horizontal_filter!(new_field, field)
    i, j, k = @index(Global, NTuple)

    new_field[i, j, k] = field[i, j, k]
    nn = (field[i, j, k], field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k])
    non_null  = Int.(nn .!= 0)

    if sum(non_null) > 0
        new_field[i, j, k] = sum(nn) / sum(non_null)
    end
end

function horizontal_filter!(new_field, field) 
    grid = field.grid
    arch = architecture(grid)

    event = launch!(arch, grid, :xyz, _horizontal_filter!, new_field, field)
    wait(Oceananigans.Architectures.device(arch), event)

    return nothing
end

@kernel function _fix_max_val!(field, max_val)
    i, j, k = @index(Global, NTuple)

    if abs(field[i, j, k]) > max_val
        field[i, j, k] = 0.0
    end
end

fix_max_val!(field, ::Nothing) = nothing

function fix_max_val!(field, max_val) 
    grid = field.grid
    arch = architecture(grid)

    event = launch!(arch, grid, :xyz, _fix_max_val!, field, max_val)
    wait(Oceananigans.Architectures.device(arch), event)

    return nothing
end

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

    event = launch!(arch, new_grid, :xyz, _horizontal_interpolate!, new_field, old_field, new_grid, old_grid, location)
    wait(Oceananigans.Architectures.device(arch), event)

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

for year in 1995:2022
    for month in 1:12
        for day in monthly_days(year)[month]
            push!(iters, parse(Int, string(year) * string(month, pad=2) * string(day, pad=2)))
        end
    end
end

# Start from day 10!!! (That is my initial condition)
iters = iters[10:end]

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
end

Nx, Ny = resolution .* (360, 150)

function read_and_interpolate_quarter_flux(name, iterations; 
                                           Nx = 360, Ny = 150,
                                           max_val = nothing, 
                                           location = (Center, Center, Center),
                                           smoothing_passes = 5)

    grid_cs510 = LatitudeLongitudeGrid(arch, size = (1440, 720, 1), latitude = (-90, 90), longitude = (-180, 180), z = (0, 1)) 
    cs_field = Field{location...}(grid_cs510)

    interp = if location[2] == Face
        zeros(Nx, Ny+1, 6)
    else
        zeros(Nx, Ny,   6)
    end

    latitude = (-75, 75) 
                            
    grid_new = LatitudeLongitudeGrid(arch; size = (Nx, Ny, 1), latitude, longitude = (-180, 180), z = (0, 1)) 
    my_field = Field{location...}(grid_new)
    my_tmp   = Field{location...}(grid_new)

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
            
        for step in 1:smoothing_passes
            horizontal_filter!(my_tmp,   my_field)
            horizontal_filter!(my_field, my_tmp)
        end

        interp[:, :, idx] .= Array(interior(my_field))
    end

    return interp
end

function transpose_flux!(var, tmp)
    Nx = size(tmp, 1)

    tmp                 .= var[Nx+1:end, :, :]
    var[Nx+1:end, :, :] .= var[1:Nx, :, :]
    var[1:Nx, :, :]     .= tmp

    return nothing
end

tmp  = zeros(Nx ÷ 2, Ny,   6)
tmpy = zeros(Nx ÷ 2, Ny+1, 6)

for (idx, iterations) in enumerate(it_collection[1:45])
    τx = read_and_interpolate_quarter_flux("oceTAUX",  iterations; smoothing_passes = 15, max_val = 1e8, location = (Face, Center, Center))
    τy = read_and_interpolate_quarter_flux("oceTAUY",  iterations; smoothing_passes = 15, max_val = 1e8, location = (Center, Face, Center))
    Fs = read_and_interpolate_quarter_flux("oceFWflx", iterations; max_val = 1e8)
    Qs = read_and_interpolate_quarter_flux("oceQnet",  iterations; max_val = 1e8)

    @info "Correcting fluxes"
    Qs .*= (1 / 1000 / 3991)
    Fs .*= (1 / 1000 * 35)

    Qs .= - Qs
    τx .= - τx ./ 1000
    τy .= - τy ./ 1000

    transpose_flux!(τx, tmp)
    transpose_flux!(Qs, tmp)
    transpose_flux!(Fs, tmp)
    transpose_flux!(τy, tmpy)

    @info "saving down fluxes $idx"
    jldsave("fluxes_$(idx).jld2", τx = τx, τy = τy, Qs = Qs, Fs = Fs)
end

