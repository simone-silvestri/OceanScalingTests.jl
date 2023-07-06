#!/bin/bash

# Let's first specify the enviromental variables
export RESOLUTION=10
export LOADBALANCE=1
export NZ=100
export EXPERIMENT="RealisticOcean"
export PROFILE=0
export WITHFLUXES=1
export FINALYEAR=1995
export FINALMONTH=2

JULIA=julia

#####
##### We need to first generate a bathymetry
#####

cd bathymetry

BATHYMETRY="bathymetry$RESOLUTION.jld2"

cat > generate_bathymetry.jl << EoF_s
include("GenerateBathymetry.jl")
using .GenerateBathymetry
res = parse(Int, get(ENV, "RESOLUTION", "3"))
bat = interpolate_bathymetry_from_ETOPO1(res, 75; interpolation_method = LinearInterpolation(passes = 5))
write_bathymetry_to_file(res, bat)
EoF_s

# check that the bathymetry exists and regenerate it if needed
if test -f "$BATHYMETRY"; then
   echo "the bathymetry file already exists."
else
    echo "regenerating bathymetry"
    JULIA --project --check-bounds=no generate_bathymetry.jl
fi

rm generate_bathymetry.jl

#####
##### Now we need to first generate daily fluxes
#####

cd ../fluxes

# check fluxes exist up to the year/month we want to simulate
FINALFLUXINDEX=$(((FINALYEAR-1995)*73+FINALMONTH*6+1))

echo $FINALFLUXINDEX

cat > write_fluxes.jl << EoF_s
include("generate_fluxes.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
generate_fluxes(parse(Int, get(ENV, "RESOLUTION", "3")); arch = CPU())
EoF_s

if test -f "fluxes_$FINALFLUXINDEX.jld2"; then
    echo "flux files already exist"
else
    echo "regenerating fluxes"
    JULIA --project --check-bounds=no write_fluxes.jl
fi

rm write_fluxes.jl

#####
##### Now we can finally run our simulation
#####

cd ../
### sbatch -N4 run_job.sh
