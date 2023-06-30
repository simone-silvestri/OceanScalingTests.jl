#!/bin/bash

# Let's first specify the enviromental variables
export RESOLUTION=12
export LOADBALANCE=1
export NZ=100
export EXPERIMENT="RealisticOcean"
export PROFILE=0
export WITHFLUXES=1
export FINALYEAR=1999

JULIA="/path/to/julia"

#####
##### We need to first generate a bathymetry
#####

cd bathymetry

cat > generate_bathymetry.jl << EoF_s
include("GenerateBathymetry.jl")
using .GenerateBathymetry
res = parse(Int, get(ENV, "RESOLUTION", "3"))
bat = interpolate_bathymetry_from_ETOPO1(res, 75; interpolation_method = LinearInterpolation(passes = 10))
write_bathymetry_to_file(bat, res)
EoF_s

JULIA --project --check-bounds=no generate_bathymetry.jl
rm generate_bathymetry.jl

#####
##### Now we need to first generate daily fluxes
#####

cd ../fluxes

cat > write_fluxes.jl << EoF_s
include("generate_fluxes.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
generate_fluxes(parse(Int, get(ENV, "RESOLUTION", "3")))
EoF_s

JULIA --project --check-bounds=no write_fluxes.jl
rm write_fluxes.jl

#####
##### Now we can finally run our simulation
#####

sbatch -N4 run_job.sh
