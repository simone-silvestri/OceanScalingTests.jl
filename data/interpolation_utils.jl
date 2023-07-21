using OceanScalingTests
using Oceananigans 
using Oceananigans.BoundaryConditions
using Oceananigans.Utils: launch!
using Oceananigans.Architectures
using Oceananigans.Architectures: architecture
using Oceananigans.Distributed
using Oceananigans.Grids: node
using Oceananigans.BoundaryConditions
using Oceananigans.BuoyancyModels: buoyancy_perturbationᶜᶜᶜ
using KernelAbstractions: @kernel, @index, @synchronize
using KernelAbstractions.Extras.LoopInfo: @unroll

@kernel function _propagate_field!(field, tmp_field, Nx, Ny)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        tmp_field[i, j, k] = field[i, j, k]

        nw = ifelse(i == 1 , field[Nx, j, k],   field[i - 1, j, k])
        ns = ifelse(j == 1 , field[i, 2, k],    field[i, j - 1, k])
        ne = ifelse(i == Nx, field[1, j, k],    field[i + 1, j, k])
        nn = ifelse(j == Ny, field[Ny-1, 2, k], field[i, j + 1, k])
        nb = (nw, ne, nn, ns)

        pw = ifelse(isnan(nw), false, true) 
        ps = ifelse(isnan(ns), false, true) 
        pe = ifelse(isnan(ne), false, true) 
        pn = ifelse(isnan(nn), false, true) 
        pb = (pw, ps, pe, pn)

        tmp_field[i, j, k] = dot(pb, nb) / sum(pb) 
    end
end

@kernel function _substitute_nans!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    if isnan(field[i, j, k]) 
        @inbounds field[i, j, k] = tmp_field[i, j, k]
    end
    if abs(field[i, j, k]) > 1e10
        @inbounds field[i, j, k] = tmp_field[i, j, k]
    end
end

function propagate_field!(field, tmp_field) 
    passes  = Ref(0)

    substitute_zeros_with_nans!(field)

    while isnan(sum(parent(field)))
        launch!(architecture(field.grid), field.grid, :xyz, _propagate_field!, field, tmp_field, field.grid.Nx, field.grid.Ny)
        launch!(architecture(field.grid), field.grid, :xyz, _substitute_nans!, field, tmp_field)
        passes[] += 1
        @info "propagate pass $(passes[]) with sum $(sum(parent(field)))"
    end

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

extend_vertically!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _extend_vertically!, field, size(field.grid, 3))

@kernel function _resort_vertically!(T, S, b, Nz, counter)
    i, j = @index(Global, NTuple)

    @inbounds begin
        counter[i, j] = 0
        @unroll for k in 2:Nz
            if b[i, j, k] < b[i, j, k-1]
                temp = T[i, j, k]
                T[i, j, k]   = T[i, j, k-1]
                T[i, j, k-1] = temp
                temp = S[i, j, k]
                S[i, j, k]   = S[i, j, k-1]
                S[i, j, k-1] = temp
                counter[i, j] += 1
            end
        end
    end
end

function resort_vertically!(T, S, buoyancy)
    grid = T.grid
    arch = architecture(grid)

    passes  = Ref(0)
    counter = arch_array(arch, ones(size(grid, 1), size(grid, 2)))

    while sum(counter) != 0
        b = KernelFunctionOperation{Center, Center, Center}(buoyancy_perturbationᶜᶜᶜ, grid, buoyancy, (; T, S))
        b = compute!(Field(b))
        launch!(architecture(grid), grid, :xy, _resort_vertically!, T, S, b, size(grid, 3), counter)
        passes[] += 1
        @info "resorting pass $(passes[])"
    end
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

horizontal_filter!(new_field, field) =
    launch!(architecture(field.grid), field.grid, :xyz, _horizontal_filter!, new_field, field)

@kernel function _cap_minimum!(field, min_val)
    i, j, k = @index(Global, NTuple)

    if field[i, j, k] < min_val
        field[i, j, k] = min_val
    end
end

cap_minimum!(field, ::Nothing) = nothing

cap_minimum!(field, min_val) =
    launch!(architecture(field.grid), field.grid, :xyz, _cap_minimum!, field, min_val)
    

@kernel function _substitute_zeros_with_nans!(field)
    i, j, k = @index(Global, NTuple)

    if field[i, j, k] == 0
        field[i, j, k] = NaN
    end
end
    
substitute_zeros_with_nans!(field) =
    launch!(architecture(field.grid), field.grid, :xyz, _substitute_zeros_with_nans!, field)

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
