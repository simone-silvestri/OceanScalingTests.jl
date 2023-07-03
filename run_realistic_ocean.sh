#!/bin/bash

# Let's first specify the enviromental variables
export RESOLUTION=12
export LOADBALANCE=1
export NZ=100
export EXPERIMENT="RealisticOcean"
export PROFILE=0
export WITHFLUXES=1
export FINALYEAR=1999
export FINALMONTH=12

JULIA="/path/to/julia"

#####
##### We need to first generate a bathymetry
#####

cd bathymetry

# check that the bathymetry exists

BATHYMETRY="bathymetry$(RESOLUTION).jld2"

if test -f "$BATHYMETRY"; then
    echo "the bathymetry file already exists."
else
    echo "regenerating bathymetry"
    cat > generate_bathymetry.jl << EoF_s
    include("GenerateBathymetry.jl")
    using .GenerateBathymetry
    res = parse(Int, get(ENV, "RESOLUTION", "3"))
    bat = interpolate_bathymetry_from_ETOPO1(res, 75; interpolation_method = LinearInterpolation(passes = 10))
    write_bathymetry_to_file(bat, res)
    EoF_s
  
    JULIA --project --check-bounds=no generate_bathymetry.jl
    rm generate_bathymetry.jl
fi
#####
##### Now we need to first generate daily fluxes
#####

cd ../fluxes

# check fluxes exist up to the year/month we want to simulate
FINALFLUXINDEX=$(((FINALYEAR-1995)*73+FINALMONTH*6+1))

if test -f "fluxes_$FINALFLUXINDEX.jld2"; then
    echo "flux files already exist."
else
    cat > write_fluxes.jl << EoF_s
    include("generate_fluxes.jl")
    res = parse(Int, get(ENV, "RESOLUTION", "3"))
    generate_fluxes(parse(Int, get(ENV, "RESOLUTION", "3")))
    EoF_s

    JULIA --project --check-bounds=no write_fluxes.jl
    rm write_fluxes.jl
fi

#####
##### Now we can finally run our simulation
#####

cd ../
sbatch -N4 run_job.sh
