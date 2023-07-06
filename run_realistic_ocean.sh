#!/bin/bash

# Let's first specify the enviromental variables

# Grid size
export RESOLUTION=12
export NZ=100

# Experiment details
export EXPERIMENT="RealisticOcean"
export WITHFLUXES=1

# Might increase performance when using a lot of cores (i.e. improves scalability)
export LOADBALANCE=1

# Do we want to profile or run a simulation?
export PROFILE=0

# How many nodes are we running on?
export NNODES=1

# The simulation always starts from 01/01/1995
export FINALYEAR=1995
export FINALMONTH=1

# Restart from interpolated fields ? "" : "numer_of_checkpointer_to_restart_from"
export RESTART=""

# Shall we regenerate fluxes?
export REGENERATEFLUXES=false

# Shall we regenerate initial conditions?
export REGENERATEINITIALCONDITIONS=true

# Julia specific enviromental variables
# export COMMON="/nobackup/users/ssilvest/perlmutter-test"
# export JULIA_DEPOT_PATH=":${COMMON}/depot"
# export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/satori"
# export JULIA_CUDA_MEMORY_POOL=none
export JULIA=julia

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
if [ test -f "$BATHYMETRY" ] || [ $EXPERIMENT -ne "RealisticOcean" ]; then
   echo "the bathymetry file already exists or we are running the DoubleDrake experiment."
else
    echo "regenerating bathymetry"
    $JULIA --project --check-bounds=no generate_bathymetry.jl
fi

rm generate_bathymetry.jl

#####
##### Now we need to generate daily fluxes and initial conditions
#####

cd ../data

# check fluxes exist up to the year/month we want to simulate
echo $REGENERATEFLUXES

cat > write_fluxes.jl << EoF_s
include("generate_fluxes.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
generate_fluxes(res; arch = CPU())
EoF_s

if [ $REGENERATEFLUXES ] || [ $WITHFLUXES -eq 0 ] || [ $EXPERIMENT -ne "RealisticOcean" ]; then
    echo "flux files already exist or we are running without fluxes"
else
    echo "regenerating fluxes"
    $JULIA --project --check-bounds=no write_fluxes.jl
fi

rm write_fluxes.jl

# check initial conditions exist otherwise regrid
cat > generate_initial_conditions.jl << EoF_s
include("regrid_initial_conditions.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
Nz  = parse(Int, get(ENV, "NZ", "100"))
regrid_initial_conditions(res, Nz; arch = CPU())
EoF_s

if [ $REGENERATEINITIALCONDITIONS ] || [ $RESTART -ne "" ] || [ $EXPERIMENT -ne "RealisticOcean" ]; then
    echo "initial condition files already exist or restarting from checkpoints"
else
    echo "regenerating initial conditions"
    $JULIA --project --check-bounds=no generate_initial_conditions.jl
fi

rm generate_initial_conditions.jl

#####
##### Now we can finally run our simulation
#####

cd ../
sbatch -N ${NNODES} satori_job.sh
