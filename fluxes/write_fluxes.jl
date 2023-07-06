include("generate_fluxes.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
generate_fluxes(parse(Int, get(ENV, "RESOLUTION", "3")); arch = CPU())
