using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.ImmersedBoundaries: immersed_cell

"""
    function load_balanced_grid(arch, N, latitude, z_faces, 
        ::Val{balance}, ::Val{experiment}) where {balance, experiment}

returns a grid which is balanced as much as possible, i.e. we try as much as possible to have a similar
number of active cells in all the ranks. This is not always possible, especially for a low number
of GPUs where the memory is already full when all Nx are the same.

The load balancing is activated id `balance == true` and `experiment == :RealisticOcean`
In any case, the computational grid has an active cells map to elide computations in the immersed domain
"""
function load_balanced_grid(arch, N, latitude, z_faces, 
                            ::Val{balance}, ::Val{experiment}) where {balance, experiment}

    Nx, Ny, Nz = N
    Nx = Nx ÷ arch.ranks[1]

    @show underlying_grid = LatitudeLongitudeGrid(arch;
                                size = (Nx, Ny, Nz),
                                longitude = (-180, 180),
                                latitude = latitude,
                                halo = (5, 5, 5),
                                z = z_faces)

    return experiment == :RealisticOcean ? 
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(realistic_bathymetry(underlying_grid)), true) :
           experiment == :DoubleDrake ?
           ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(double_drake_bathymetry)) :
           underlying_grid
end

function load_balanced_grid(arch, N, latitude, z_faces, 
                            ::Val{true}, ::Val{:RealisticOcean}) 

    underlying_grid = LatitudeLongitudeGrid(GPU();
                        size = N,
                        longitude = (-180, 180),
                        latitude = latitude,
                        halo = (5, 5, 5),
                        z = z_faces)

    ibg = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(realistic_bathymetry(underlying_grid))) 

    load_per_slab = arch_array(GPU(), zeros(Int, N[1]))

    loop! = assess_load(device(GPU()), 512, N[1])
    loop!(load_per_slab, ibg)

    load_per_slab = arch_array(CPU(), load_per_slab)
    local_N = calculate_local_N(load_per_slab, N, arch)

    # We cannot have Nx > 650 if Nranks = 32 otherwise we incur in memory limitations,
    # so for a small number of GPUs we are limited in the load balancing
    redistribute_size_to_fulfill_memory_limitation!(local_N, 650)

    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    N    = (local_N[rank+1], N[2], N[3])

    @info "slab decomposition with " rank N

    @show underlying_grid = LatitudeLongitudeGrid(arch;
                                                  size = N,
                                                  longitude = (-180, 180),
                                                  latitude = latitude,
                                                  halo = (5, 5, 5),
                                                  z = z_faces)

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(realistic_bathymetry(underlying_grid)), true) 
end

@kernel function assess_load(load_per_slab, ibg)
    i = @index(Global, Linear)

    @unroll for j in 1:size(ibg, 2)
        @unroll for k in 1:size(ibg, 3)
            @inbounds load_per_slab[i] += ifelse(immersed_cell(i, j, k, ibg), 0, 1)
        end
    end
end

function calculate_local_N(load_per_slab, N, arch)
    active_cells  = sum(load_per_slab)
    active_load   = active_cells / arch.ranks[1]
    local_N = zeros(Int, arch.ranks[1]) # fill the local N with the active load
    idx = 1
    for r in 1:arch.ranks[1]-1
        local_load = 0
        while local_load <= active_load
            local_load += load_per_slab[idx]
            local_N[r] += 1
            idx += 1
        end
    end

    local_N[end] = N[1] - sum(local_N[1:end-1])

    return local_N
end

function redistribute_size_to_fulfill_memory_limitation!(l, m = 700)
    n = length(l)
    while any(l .> m)
        x⁺, i⁺ = findmax(l)
        x⁻, i⁻ = findmin(l)
        if x⁺ == m + 1
            l[i⁺] -= 1
            l[i⁻] += 1
        else
            Δ = l[i⁺] - m
            q = Δ ÷ n
            l .+= q
            l[i⁺] -= Δ
            l[i⁻] += mod(Δ, n)
        end
    end
    while any(l .< 20)
        x⁺, i⁺ = findmax(l)
        x⁻, i⁻ = findmin(l)
        if x⁻ == 19
            l[i⁻] += 1
            l[i⁺] -= 1
        else
            Δ = 20 - l[i⁻]
            q = Δ ÷ n
            l .-= q
            l[i⁺] -= mod(Δ, n)
            l[i⁻] += Δ
        end
    end
end