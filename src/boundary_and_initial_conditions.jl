using Oceananigans.Units
using Oceananigans.Utils
using Oceananigans.Grids: halo_size
using Oceananigans.Architectures: arch_array, architecture, device
using Oceananigans.DistributedComputations: partition_global_array
using OffsetArrays
using KernelAbstractions: @index, @kernel

initialize_model!(model, ::Val{:Quiescent}; kw...)   = nothing

function initialize_model!(model, ::Val{:DoubleDrake}; kw...)
    Depth = model.grid.Lz
    exp_length = Depth / 15.0

    @inline init_temperature(λ, φ, z) = exponential_profile(z; Lz = Depth, h = exp_length) * max(T_reference(φ), 2.0)

    set!(model, T = init_temperature, S = 35.0)

    return nothing
end

function initialize_model!(model, ::Val{:RealisticOcean}; restart = "")
    
    if !isempty(restart)
        return nothing
    end

    grid = model.grid

    arch = architecture(grid)
    T_init = zeros(size(grid)...)
    S_init = zeros(size(grid)...)

    rx = grid.architecture.local_rank

    for k in 1:size(grid, 3)
       if rx == 1
          @info "loading level $k"
       end

       T_init[:, :, k] .= Array(partition_global_array(arch, jldopen("data/initial_T_at_k$(k).jld2")["T"][:, :, 1], size(grid)[[1, 2]]))
       S_init[:, :, k] .= Array(partition_global_array(arch, jldopen("data/initial_S_at_k$(k).jld2")["S"][:, :, 1], size(grid)[[1, 2]]))
    end

    set!(model, T = T_init, S = S_init)

    return nothing
end

function set_boundary_conditions(::Val{:Quiescent}, grid; kw...)

    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function set_boundary_conditions(::Val{:DoubleDrake}, grid; kw...)

    Ly = grid.Ly

    S_reference = eltype(grid)(35.0) # reference salinity for salinity flux (psu)
    u_coeffs = wind_stress_coefficients(Ly/2)
    S_coeffs = salinity_flux_coefficients(Ly/2, S_reference)

    λ = eltype(grid)(5meters / 7days) # Pumping velocity of temperature restoring (ms⁻¹)
    u_top_bc = FluxBoundaryCondition(surface_stress_x,      discrete_form=true, parameters=u_coeffs)
    S_top_bc = FluxBoundaryCondition(surface_salinity_flux, discrete_form=true, parameters=S_coeffs)
    T_top_bc = FluxBoundaryCondition(T_relaxation,          discrete_form=true, parameters=λ)

    μ = eltype(grid)(0.003) # linear drag coefficient (ms⁻¹)
    u_bot_bc = FluxBoundaryCondition(u_linear_bottom_drag, discrete_form=true, parameters=μ)
    v_bot_bc = FluxBoundaryCondition(v_linear_bottom_drag, discrete_form=true, parameters=μ)

    T_bcs = FieldBoundaryConditions(top=T_top_bc)
    S_bcs = FieldBoundaryConditions(top=S_top_bc)
    u_bcs = FieldBoundaryConditions(bottom=u_bot_bc, top=u_top_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bot_bc)

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function set_boundary_conditions(::Val{:RealisticOcean}, grid; with_fluxes = true, with_restoring = true)

    arch = architecture(grid)

    Nx, Ny, _ = size(grid)

    if with_fluxes
        τx = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny,   6), (-1, 0, 0)))
        τy = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny+1, 6), (-1, 0, 0)))
        Qs = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny,   6), (-1, 0, 0)))
        Fs = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny,   6), (-1, 0, 0)))
        
        load_fluxes!(grid, τx, τy, Qs, Fs, 1)
        
        u_top_bc = FluxBoundaryCondition(flux_from_interpolated_array, discrete_form=true, parameters=τx)    
        v_top_bc = FluxBoundaryCondition(flux_from_interpolated_array, discrete_form=true, parameters=τy)
        if with_restoring
            Tr = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny, 6), (-1, 0, 0)))
            Sr = arch_array(arch, OffsetArray(zeros(eltype(grid), Nx+2, Ny, 6), (-1, 0, 0)))

            load_restoring!(grid, Tr, Sr, 1)

            Δz = minimum_zspacing(grid)

            # If temperature or salinity are outside the physical range, crank up restoring velocity
            T_top_bc = FluxBoundaryCondition(flux_and_restoring_T, discrete_form=true, parameters=(; Qs, Tr, λ=Δz/30days))
            S_top_bc = FluxBoundaryCondition(flux_and_restoring_S, discrete_form=true, parameters=(; Fs, Sr, λ=Δz/60days))
        else
            T_top_bc = FluxBoundaryCondition(flux_from_interpolated_array, discrete_form=true, parameters=Qs)
            S_top_bc = FluxBoundaryCondition(flux_from_interpolated_array, discrete_form=true, parameters=Fs)
        end
    else
        u_top_bc = FluxBoundaryCondition(zero(grid))    
        v_top_bc = FluxBoundaryCondition(zero(grid))
        T_top_bc = FluxBoundaryCondition(zero(grid))
        S_top_bc = FluxBoundaryCondition(zero(grid))
    end

    μ = eltype(grid)(0.001) # Quadratic drag coefficient (ms⁻¹)
    u_bot_bc = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=μ)
    v_bot_bc = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_bot_bc = FluxBoundaryCondition(u_immersed_quadratic_bottom_drag, discrete_form=true, parameters=μ)
    v_immersed_bot_bc = FluxBoundaryCondition(v_immersed_quadratic_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = u_immersed_bot_bc)
    v_immersed_bc = ImmersedBoundaryCondition(bottom = v_immersed_bot_bc)

    T_bcs = FieldBoundaryConditions(top=T_top_bc)
    S_bcs = FieldBoundaryConditions(top=S_top_bc)
    u_bcs = FieldBoundaryConditions(bottom=u_bot_bc, top=u_top_bc, immersed = u_immersed_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bot_bc, top=v_top_bc, immersed = v_immersed_bc)

    return (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs)
end

function load_fluxes!(grid, τx, τy, Qs, Fs, filenum)
    arch = architecture(grid)
    rx   = arch.local_rank
    nx, ny, _ = size(grid) 

    if rx == 1
        @info "loading fluxes from file fluxes_$(filenum).jld2"
    end

    file = jldopen("data/fluxes_$(filenum).jld2")
    
    # Garbage collect!!
    GC.gc()

    τxtmp = partition_global_array(arch, file["τx"], (nx, ny,   6))
    τytmp = partition_global_array(arch, file["τy"], (nx, ny+1, 6))
    Qstmp = partition_global_array(arch, file["Qs"], (nx, ny,   6))
    Fstmp = partition_global_array(arch, file["Fs"], (nx, ny,   6))

    _copy_flux_values(device(arch), (16, 16), (nx, ny, 6))(τx, τy, Qs, Fs, τxtmp, τytmp, Qstmp, Fstmp, nx, ny)

    return nothing
end

@kernel function _copy_flux_values!(τx, τy, Qs, Fs, τxtmp, τytmp, Qstmp, Fstmp, Nx, Ny)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        τx[i, j, k] = τxtmp[i, j, k]
        τy[i, j, k] = τytmp[i, j, k]
        Qs[i, j, k] = Qstmp[i, j, k]
        Fs[i, j, k] = Fstmp[i, j, k]
    end

    if j == 1
        τy[i, Ny+1, k] = τytmp[i, Ny+1, k]
    end

    if i == 1
        τx[0, j, k]    = τxtmp[Nx, j, k]
        τx[Nx+1, j, k] = τxtmp[1, j, k]
        τy[0, j, k]    = τytmp[Nx, j, k]
        τy[Nx+1, j, k] = τytmp[1, j, k]
        Qs[0, j, k]    = Qstmp[Nx, j, k]
        Qs[Nx+1, j, k] = Qstmp[1, j, k]
        Fs[0, j, k]    = Fstmp[Nx, j, k]
        Fs[Nx+1, j, k] = Fstmp[1, j, k]
    end
end

function load_restoring!(grid, Tr, Sr, filenum)
    arch = architecture(grid)
    rx   = arch.local_rank
    nx, ny, _ = size(grid) 

    if rx == 1
        @info "loading fluxes from file restoring_$(filenum).jld2"
    end

    file = jldopen("data/restoring_$(filenum).jld2")
    
    # Garbage collect!!
    GC.gc()

    Trtmp = partition_global_array(arch, file["Tr"], (nx, ny, 6))
    Srtmp = partition_global_array(arch, file["Sr"], (nx, ny, 6))

    _copy_restoring_values(device(arch), (16, 16), (nx, ny, 6))(Tr, Sr, Trtmp, Srtmp, nx, ny)

    return nothing
end

@kernel function _copy_restoring_values(Tr, Sr, Trtmp, Srtmp, nx, ny)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Tr[i, j, k] = Trtmp[i, j, k]
        Sr[i, j, k] = Srtmp[i, j, k]
    end

    if i == 1
        Tr[0, j, k]    = Trtmp[Nx, j, k]
        Tr[Nx+1, j, k] = Trtmp[1, j, k]
        Sr[0, j, k]    = Srtmp[Nx, j, k]
        Sr[Nx+1, j, k] = Srtmp[1, j, k]
    end
end

# Update fluxes through a Callback every 5days.
@inline function update_fluxes!(sim)

    # Repeat year does mod(time, 365) otherwise take out the mod
    repeat_year = parse(Bool, get(ENV, "REPEATYEAR", "true"))
    filenum = repeat_year ? Int(mod(sim.model.clock.time ÷ days, 365) ÷ 5) + 1 :
                                Int(sim.model.clock.time ÷ 5days) + 1 

    model = sim.model
    grid  = model.grid

    τx = model.velocities.u.boundary_conditions.top.condition.parameters 
    τy = model.velocities.v.boundary_conditions.top.condition.parameters 

    Qs = model.tracers.T.boundary_conditions.top.condition.parameters 
    Fs = model.tracers.S.boundary_conditions.top.condition.parameters 

    # Take care of case with restoring
    Qs = Qs isa NamedTuple ? Qs.Qs : Qs
    Fs = Fs isa NamedTuple ? Fs.Fs : Fs

    load_fluxes!(grid, τx, τy, Qs, Fs, filenum)
    
    return nothing
end

# Update restoring through a Callback every 15days.
@inline function update_restoring!(sim)

    # Repeat year does mod(time, 365) otherwise take out the mod
    repeat_year = parse(Bool, get(ENV, "REPEATYEAR", "true"))
    if repeat_year
        year_days = mod(sim.model.clock.time ÷ days, 365)
        filenum = mod(year_days, 15) == 0 ? Int(year_days ÷ 15) + 1 : nothing
    else
        filenum = Int(sim.model.clock.time ÷ 15days) + 1 
    end

    if !isnothing(filenum)
        model = sim.model
        grid  = model.grid

        Tr = model.tracers.T.boundary_conditions.top.condition.parameters.Tr 
        Sr = model.tracers.S.boundary_conditions.top.condition.parameters.Sr 

        load_restoring!(grid, Tr, Sr, filenum)
    end
    
    return nothing
end
