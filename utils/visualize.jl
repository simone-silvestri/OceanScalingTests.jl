using JLD2
using GLMakie

function geographic2cartesian(λ, φ; r=1)
    Nλ = length(λ)
    Nφ = length(φ)

    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ) 
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

function visualize_makie(vars, λs, ϕ, iter, iters, scale, map)

    fig = Figure(resolution = (700, 700))
    ax = fig[:, :] = LScene(fig, show_axis=false) # make plot area wider

    for (var, λ) in zip(vars, λs)
        x, y, z = geographic2cartesian(λ, ϕ, r=1.01)
        surface!(ax, x, y, z, color=var, colorrange=scale, colormap=map) 
    end

    init = (π/5, π, 0)
    rotate_cam!(ax.scene, init)
    rot  = (0, π/100, 0)

    record(fig, "output.mp4", iters, framerate=10) do i
        @info "Plotting iteration $i of $(iters[end])..."
        rotate_cam!(ax.scene, rot)
        iter[] = i
    end
    return nothing
end

function visualize_globe(var, prefix; 
                         Nranks = 2, 
                         resolution = 3,
                         colorrange = (-0.5e-5, 0.5e-5), 
                         colormap = :solar,
                         halos = (5, 5))

    Nx = 360 * resolution
    nx = Nx ÷ Nranks
    Hx, Hy = halos

    files = [jldopen("$(prefix)_$(rank).jld2")   for rank in 0:Nranks-1]
    λs    = [(1+rank*nx:(rank+1)*nx) .* 360 / Nx for rank in 0:Nranks-1]
    
    φ = collect(-75+1/resolution:1/resolution:75)

    iter = Observable(0)

    vs = [@lift(files[r]["timeseries/$(var)/" * string($iter)][Hx+1:end-Hx, Hy+1:end-Hy, 1]) for r in 1:Nranks]

    iters = parse.(Int, keys(files[1]["timeseries/t"]))

    visualize_makie(vs, λs, φ, iter, iters, colorrange, colormap)

    for rank in 1:Nranks
        close(files[rank])
    end
end
