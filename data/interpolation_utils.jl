using OceanScalingTests
using Oceananigans 
using Oceananigans.BoundaryConditions
using Oceananigans.Utils: launch!
using Oceananigans.Architectures
using Oceananigans.Architectures: architecture
using Oceananigans.Distributed
using Oceananigans.Grids: node
using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

@kernel function _propagate_field!(field, Nx, Ny)
    i, j, k = @index(Global, NTuple)
    nw = ifelse(i == 1 , field[Nx, j, k], field[i - 1, j, k])
    ns = ifelse(j == 1 , 0.0,             field[i, j - 1, k])
    ne = ifelse(i == Nx, field[1, j, k],  field[i + 1, j, k])
    nn = ifelse(j == Ny, 0.0,             field[i, j + 1, k])
    nb = (nw, ne, nn, ns)
    pos = Int.(nb .!= 0)

    @synchronize 

    if (field[i, j, k] == 0) & (sum(pos) > 0)
        field[i, j, k] = dot(pos, nb) / sum(pos) 
    end
end

propagate_field!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _propagate_field!, field, field.grid.Nx, field.grid.Ny)

@kernel function _extend_vertically!(field, Nz)
    i, j = @index(Global, NTuple)

    @unroll for k in Nz-1:-1:1
        if field[i, j, k] == 0.0
            field[i, j, k] = field[i, j, k+1]
        end
    end
end

extend_vertically!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _extend_vertically!, field, size(field.grid, 3))

@kernel function _horizontal_filter!(new_field, field)
    i, j, k = @index(Global, NTuple)

    new_field[i, j, k] = field[i, j, k]
    nn = (field[i, j, k], field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k])
    non_null  = Int.(nn .!= 0)

    @synchronize

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
    x, y, z = node(i, j, k, new_grid, loc...)
    
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
