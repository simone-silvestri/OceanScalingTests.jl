using Oceananigans.Units
using Oceananigans.Grids: topology
using Oceananigans.Operators: assumed_field_location
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

function check_ranges(folder, ranks; H = 7, fields = false, time = 0, iteration = 0)
    Nx = Vector(undef, ranks)
    iranges = Vector(undef, ranks)
    for rank in 0:ranks - 1
        var = fields ? jldopen(folder * "RealisticOcean_fields_$(rank).jld2")["timeseries/u/"*string(time)][H+1:end-H, H+1:end-H, 1] :
                       jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["u/data"][H+1:end-H, H+1:end-H, H+1:end-H] 
        Nx[rank+1] = size(var, 1)
    end
    iranges[1] = UnitRange(1, Nx[1])
    for rank in 2:ranks
        iranges[rank] = UnitRange(iranges[rank-1][end]+1,iranges[rank-1][end]+Nx[rank])
    end

    return iranges
end

function compress_restart_file(full_size, ranks, iteration, folder = "../"; Depth = 5244.5, 
                               bathymetry = jldopen("data/bathymetry.jld2")["bathymetry"], H = 7)

    Nx, Ny, Nz = full_size

    z_faces = OceanScalingTests.exponential_z_faces(Nz, Depth)

    full_grid = LatitudeLongitudeGrid(; size = full_size,
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (5, 5, 5),
                                        z = z_faces)

    @info "initializing active map"
    fields_data = Dict()
    fields_data[:underlying_grid] = full_grid
    fields_data[:bathymetry]      = bathymetry
    fields_data[:clock] = jldopen(folder * "RealisticOcean_checkpoint_0_iteration$(iteration).jld2")["clock"]

    iranges = check_ranges(folder, ranks; H, iteration)

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

    compressed_η = zeros(Float32, Nx, Ny, 1)
    for rank in 0:ranks-1
        @info "reading rank $rank"

        irange = iranges[rank+1]
        data = jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["η/data"]
        Hx = calc_free_surface_halo(irange, data)
        data = data[Hx+1:end-Hx, H+1:end-H, :]
        compressed_η[irange, :, :] .= Float32.(data)
    end

    fields_data[:η] = compressed_η

    jldopen(folder * "compressed_iteration_$(iteration).jld2","w") do f
        for (key, value) in fields_data
            f[string(key)] = value
        end
    end
end

function calc_free_surface_halo(irange, data)
    Nx = size(data, 1)
    nx = length(irange)
    return Int((Nx - nx) ÷ 2)
end

function compress_surface_fields(full_size, ranks, folder = "../"; suffix = "", H = 7)

    Nx, Ny, Nz = full_size
    Nz = 1

    fields_data = Dict()

    times   = keys(jldopen(folder * "RealisticOcean_fields_0.jld2")["timeseries/t"])
    iranges = check_ranges(folder, ranks; H, fields = true, time = times[1])
    Nt = length(times)

    @info "starting the compression"
    for var in (:u, :w, :v, :T, :S)
        GC.gc()

        @info "compressing variable $var"
        sizefield = var == :v ? (Nx, Ny+1, Nt) :
                    var == :w ? (Nx, Ny  , Nt) : (Nx, Ny, Nt)

        compressed_data = zeros(Float32, sizefield...)

        for rank in 0:ranks-1
            @info "reading rank $rank"
            irange = iranges[rank+1]
            for (i, t) in enumerate(times)
                compressed_data[irange, :, i] .= jldopen(folder * "RealisticOcean_fields_$(rank).jld2")["timeseries/"*string(var)*"/"*t][H+1:end-H, H+1:end-H, 1]
            end
        end

        fields_data[var] = compressed_data
    end

    jldopen(folder * "compressed_surface_fields" * suffix * ".jld2","w") do f
        for (key, value) in fields_data
            f[string(key)] = value
        end
    end
end

const regex = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";

function compress_all_restarts(full_size, ranks, dir)
    files = readdir(dir)
    files = filter(x -> length(x) > 30, files)
    files = filter(x -> x[1:26] == "RealisticOcean_checkpoint_", files)
    iterations = Int[]
    for file in files
        file   = file[1:end-5]
        string = ""
        i = length(file)
        while occursin(regex, "$(file[i])")
            string = file[i] * string
            i -= 1
        end
        push!(iterations, parse(Int, string))
    end

    iterations = unique(iterations)
    iterations = sort(iterations)
    for iter in iterations
        @info "compressing iteration $iter"
        compress_restart_file(full_size, ranks, iter, dir; bathymetry = nothing)
    end
end


