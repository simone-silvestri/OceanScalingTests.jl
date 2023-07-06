using Oceananigans 
using Oceananigans.BoundaryConditions
using Oceananigans.Architectures: device, arch_array
using Oceananigans.Grids: architecture
using Oceananigans.Fields: regrid!, regrid_in_y!, regrid_in_x!
using Oceananigans.Utils: launch!
using DataDeps
using JLD2 
using Statistics: dot
using KernelAbstractions: @index, @kernel, @synchronize
using KernelAbstractions.Extras.LoopInfo: @unroll

include("interpolation_utils.jl")

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_z_faces(Nz, Depth; h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - Depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
end

datadep"quarter_degree_near_global_lat_lon"

datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/z_faces-50-levels.jld2"
file_z_faces = jldopen(datadep_path)
ECCO_z_faces = file_z_faces["z_faces"][3:end];

function regrid_initial_conditions(resolution, Nz; arch = GPU(), 
                                   regrid_in_z = true, regrid_in_x = true, z_faces = ECCO_z_faces,
                                   filepath = "evolved-initial-conditions-1105days.jld2")

    file_init = jldopen(filepath)

    # old grid size 
    nx, ny, nz = size(file_init["T"])

    # new grid size
    Nx = Int(360 * resolution) 
    Ny = Int(150 * resolution)

    # the initial grid
    grid = LatitudeLongitudeGrid(size = (nx, ny, nz),
                                 longitude = (-180, 180),
                                 latitude = (-75, 75),
                                 z = z_faces)

    T = set!(CenterField(grid), file_init["T"])
    S = set!(CenterField(grid), file_init["S"])

    @info "extending vertically T and S"
    extend_vertically!(T)
    extend_vertically!(S)

    @info "finished extending"

    Depth = grid.Lz

    @info "start regridding in Z!!"
    gridᶻ = LatitudeLongitudeGrid(size = (nx, ny, Nz),
                                  longitude = (-180, 180),
                                  latitude = (-75, 75),
                                  z = exponential_z_faces(Nz, Depth))

    Tᶻ = CenterField(gridᶻ)
    Sᶻ = CenterField(gridᶻ)

    if regrid_in_z
        fill_halo_regions!((T, S))
        regrid!(Tᶻ, T)
        regrid!(Sᶻ, S)

        for step in 1:50
            @info "propagating step $step"
            propagate_field!(Tᶻ)
            propagate_field!(Sᶻ)
        end

        jldsave("regridded_in_z.jld2", T = Array(interior(Tᶻ)), S = Array(interior(Sᶻ)))
    else
        set!(Tᶻ, jldopen("regridded_in_z.jld2")["T"])
        set!(Sᶻ, jldopen("regridded_in_z.jld2")["S"])
    end

    @info "Finished regridding in z"

    # The regridding will now be done per level to avoid OOM errors
    gridᶻ¹ = LatitudeLongitudeGrid(arch, size = (nx, ny, 1),
                                  longitude = (-180, 180),
                                  latitude = (-75, 75),
                                  z = (0, 1))
    
    Tᶻ¹ = CenterField(gridᶻ¹)
    Sᶻ¹ = CenterField(gridᶻ¹)
                                 
    @info "Continue by regridding in X!!"
    gridˣᶻ = LatitudeLongitudeGrid(arch, size = (Nx, ny, 1),
                                   longitude = (-180, 180),
                                   latitude = (-75, 75),
                                   z = (0, 1))

    Tˣᶻ = CenterField(gridˣᶻ)
    Sˣᶻ = CenterField(gridˣᶻ)

    if regrid_in_x
        for k in Nz:-1:1
            @info "regridding level $(k)"

            set!(Tᶻ¹, arch_array(arch, interior(Tᶻ, :, :, k:k)))
            set!(Sᶻ¹, arch_array(arch, interior(Sᶻ, :, :, k:k)))

            fill_halo_regions!((Tᶻ¹, Sᶻ¹))

            regrid!(Tˣᶻ, Tᶻ¹)
            regrid!(Sˣᶻ, Sᶻ¹)

            jldsave("T_regridded_in_xz_at$(k).jld2", T = Array(interior(Tˣᶻ)))
            jldsave("S_regridded_in_xz_at$(k).jld2", S = Array(interior(Sˣᶻ)))
        else
    end

    @info "Conclude by regridding in Y!!"
    gridˣʸᶻ = LatitudeLongitudeGrid(arch, size = (Nx, Ny, 1),
                                    longitude = (-180, 180),
                                    latitude = (-75, 75),
                                    z = (0, 1))
            
    Tˣʸᶻ = CenterField(gridˣʸᶻ)
    Sˣʸᶻ = CenterField(gridˣʸᶻ)                 

    for k in 1:Nz        
        @info "regridding level $(k)"

        set!(Tˣᶻ, jldopen("T_regridded_in_xz_at$(k).jld2")["T"])
        set!(Sˣᶻ, jldopen("S_regridded_in_xz_at$(k).jld2")["T"])

        fill_halo_regions!((Tˣᶻ, Sˣᶻ))
        
        regrid!(Tˣʸᶻ, Tˣᶻ)
        regrid!(Sˣʸᶻ, Sˣᶻ)

        jldsave("initial_T_at_k$(k).jld2", T = Array(interior(Tˣʸᶻ)))
        jldsave("initial_S_at_k$(k).jld2", S = Array(interior(Sˣʸᶻ)))
    end
end