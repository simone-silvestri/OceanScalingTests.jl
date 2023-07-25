using Oceananigans.Units
using Oceananigans.Grids: topology
using Oceananigans.Operators: assumed_field_location
using Oceananigans.ImmersedBoundaries: ActiveCellsIBG
using KernelAbstractions: @kernel, @index
using JLD2

# Doing the year 1995-2020
leap_year_days(year) = year == 1996 || 
                       year == 2000 || 
                       year == 2004 || 
                       year == 2008 || 
                       year == 2012 || 
                       year == 2016 || 
                       year == 2020 ? UnitRange(1, 29) : UnitRange(1, 28)

monthly_days(year) = [1:31, leap_year_days(year), 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

function realistic_ocean_stop_time(final_year, final_month = 12)
    simulation_days = 0
    # Starts from 01/01/1995
    for year in 1995:final_year-1
        days_in_month = monthly_days(year)
	simulation_days += sum(days_in_month[m][end] for m in 1:12)
    end

    days_in_month = monthly_days(final_year)
    simulation_days += sum(days_in_month[m][end] for m in 1:final_month)

    return simulation_days * days 
end

function check_ranges(folder, ranks, iteration; H = 7)
    Nx = Vector(undef, ranks)
    iranges = Vector(undef, ranks)
    for rank in 0:ranks - 1
        var = jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["u/data"][H+1:end-H, H+1:end-H, H+1:end-H]
        Nx[rank+1] = size(var, 1)
    end
    iranges[1] = UnitRange(1, Nx[1])
    for rank in 2:ranks
        iranges[rank] = UnitRange(iranges[rank-1][end]+1,iranges[rank-1][end]+Nx[rank])
    end

    return iranges
end

function compress_restart_file(full_size, ranks, iteration, folder = "../"; Depth = 5244.5, 
                               bathymetry = jldopen("data/bathymetry.jld2")["bathymetry"], Nsteps = 33, H = 7)

    Nx, Ny, Nz = full_size

    z_faces = OceanScalingTests.exponential_z_faces(Nz, Depth)

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

    iranges = check_ranges(folder, ranks, iteration; H)

    @info "starting the compression"
    for var in (:u, :w, :v, :T, :S)
        GC.gc()

        @info "compressing variable $var"
        sizefield = var == :v ? (Nx, Ny+1, Nz) :
                    var == :w ? (Nx, Ny, Nz+1) : (Nx, Ny, Nz)

        compressed_data = zeros(Float32, sizefield...)

        for rank in 0:ranks-1
            @info "reading rank $rank"
            irange = iranges[rank+1]
            compressed_data[irange, :, :] .= jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")[string(var) * "/data"][H+1:end-H, H+1:end-H, H+1:end-H]
        end

        fields_data[var] = compressed_data
    end

    # compressed_η = zeros(Float32, Nx, Ny, 1)
    # for rank in 0:ranks-1
    #     @info "reading rank $rank"

    #     irange = iranges[rank+1]
    #     data = jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["η/data"][Nsteps+1:end-Nsteps, 6:end-5, :]
    #     compressed_η[irange, :, :] .= Float32.(data)
    # end

    # fields_data[:η] = compressed_η

    jldopen("compressed_iteration_$(iteration).jld2","w") do f
        for (key, value) in fields_data
            f[string(key)] = value
        end
    end
end