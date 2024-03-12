using JLD2
using CairoMakie

const regex = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";

files = readdir("./")
files = filter(x -> length(x) > 20, files)
files = filter(x -> x[1:20] == "compressed_iteration", files)
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

iter = Observable(iterations[1])

bath = jldopen("../bathymetry/bathymetry12.jld2")["bathymetry"]

e = @lift begin
     file = jldopen("compressed_iteration_" * string($iter) * ".jld2")
     e1 = file["e"][:, :, 100]
     e1[bath .> 0] .= NaN;
     e1
end

fig = Figure(resolution = (1000, 500))

ax = Axis(fig[1, 1])
heatmap!(e, colorrange = (0, 0.001))
hidedecorations!(ax)
hidespines!(ax)

j = Ref(1)

CairoMakie.record(fig, "tke.mp4", iterations) do i 
	iter[] = i
	@info "doing iteration $i"
	CairoMakie.save("iteration_new$(j[]).png", fig)
	j[] = j[] + 1
end
