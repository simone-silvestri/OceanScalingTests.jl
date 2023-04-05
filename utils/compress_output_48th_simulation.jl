using Oceananigans
using Oceananigans.Grids: topology
using Oceananigans.Operators: assumed_field_location
using Oceananigans.ImmersedBoundaries: ActiveCellsIBG, active_cells_map
using Oceananigans.Grids: architecture
using Oceananigans.Architectures: device
using OceanScalingTests
using BitInformation: round!
using KernelAbstractions: @kernel, @index
using JLD2

function compress_restart_file(full_size, ranks, iteration, folder = "../"; mantissa_bits = 15, Depth = 5244.5, bathymetry = jldopen("data/bathymetry.jld2")["bathymetry"])

    Nx, Ny, Nz = full_size

    z_faces = OceanScalingTests.exponential_z_faces(Nz, Depth)
    nx = Nx ÷ ranks

    @show size

    full_grid = LatitudeLongitudeGrid(; size = full_size,
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (5, 5, 5),
                                        z = z_faces)

    @info "initializing active map"
    fields_data = Dict()
    fields_data[:underlying_grid] = full_grid
    fields_data[:bathymetry]      = bathymetry

    @info "starting the compression"
    for var in (:u, :v, :w, :T, :S)
        GC.gc()

        @info "compressing variable $var"
        sizefield = var == :v ? (Nx, Ny+1, Nz) :
                    var == :w ? (Nx, Ny, Nz+1) : (Nx, Ny, Nz)

	    compressed_data = zeros(Float32, sizefield...)

        for rank in 0:ranks-1
            @info "reading rank $rank"
            irange    = UnitRange(1 + rank * nx, (rank + 1) * nx)               
            compressed_data[irange, :, :] .= jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")[string(var) * "/data"][6:end-5, 6:end-5, 6:end-5]
        end

        # compressed_data  = compress_field_interior(compressed_data, full_grid; mantissa_bits)
        fields_data[var] = compressed_data
    end

    compressed_η = zeros(Float32, Nx, Ny, 1)
    for rank in 0:ranks-1
        @info "reading rank $rank"
        
        irange = UnitRange(1 + rank * nx, (rank + 1) * nx)
        data = jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["η/data"][34:end-33, 6:end-5, :]
        compressed_η[irange, :, :] .= Float32.(data)
    end

    fields_data[:η] = compressed_η

    jldopen("compressed_iteration_$(iteration).jld2","w") do f
        for (key, value) in fields_data
            f[string(key)] = value
        end
    end
end

@kernel function fill_compressed_data!(compressed_data, field, index_list)
    idx = @index(Global, Linear)
    
    @inbounds begin
        i, j, k = index_list[idx]
        compressed_data[idx] = field[i, j, k]
    end
end

@kernel function read_compressed_data!(field, compressed_data, index_list)
    idx = @index(Global, Linear)
    
    @inbounds begin
        i, j, k = index_list[idx]
        field[i, j, k] = compressed_data[idx]
    end
end

function compress_field_interior(field::Field, grid; mantissa_bits = 23) 
    compressed_data = interior(field) |> Array{Float32}
    round!(compressed_data, mantissa_bits)

    return compressed_data
end

function compress_field_interior(data::AbstractArray{T, 3}, grid; mantissa_bits = 23) where T
    compressed_data = data |> Array{Float32}
    round!(compressed_data, mantissa_bits)

    return compressed_data
end

function compress_field_interior(field::Field, grid::ActiveCellsIBG; mantissa_bits = 23) 
    arch            = architecture(grid)
    index_list      = grid.active_cells_map
    compressed_data = zeros(Float32, length(index_list)) 

    loop! = fill_compressed_data!(device(arch), 256, length(index_list))
    loop!(compressed_data, field, index_list)

    round!(compressed_data, mantissa_bits)

    return compressed_data
end

function compress_field_interior(data::AbstractArray{T, 3}, grid::ActiveCellsIBG; mantissa_bits = 23) where T 
    arch            = architecture(grid)
    index_list      = grid.active_cells_map
    compressed_data = zeros(Float32, length(index_list)) 

    loop! = fill_compressed_data!(device(arch), 256, length(index_list))
    loop!(compressed_data, data, index_list)

    round!(compressed_data, mantissa_bits)

    return compressed_data
end

function decompress_data_into_field!(field, compressed_data, index_list)
    arch = architecture(field)
    loop! = read_compressed_data!(device(arch), 256, length(index_list))
    loop!(field, compressed_data, index_list)

    return nothing
end

N = (48 * 360, 48 * 150, 100)

compress_restart_file(N, 32, 334290, "fields-150-to-190/"; mantissa_bits = 16, bathymetry = nothing)
compress_restart_file(N, 32, 355890, "fields-150-to-190/"; mantissa_bits = 16, bathymetry = nothing)
compress_restart_file(N, 32, 377490, "fields-150-to-190/"; mantissa_bits = 16, bathymetry = nothing)
