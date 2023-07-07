#!/bin/bash

# Let's first specify the enviromental variables

# ====================================================== #
# ================ USER SPECIFIED INPUTS =============== #
# ====================================================== #

# Grid size
export RESOLUTION=12
export NZ=100

# Experiment details
export EXPERIMENT="RealisticOcean"
export WITHFLUXES=0
export PRECISION="Float64"

# Might increase performance when using a lot of cores (i.e. improves scalability)
export LOADBALANCE=1

# Do we want to profile or run a simulation?
export PROFILE=0

# How many nodes are we running on?
export NNODES=1

# For :RealisticOcean the simulation always starts from 01/01/1995, when does it end?
export FINALYEAR=1995
export FINALMONTH=1

# Restart from interpolated fields ? "" : "numer_of_iteration_to_restart_from"
export RESTART=""

# Shall we regenerate fluxes? Shall we regenerate initial conditions?
export REGENERATEFLUXES=0
export REGENERATEINITIALCONDITIONS=1

# Julia specific enviromental variables
# export COMMON="/nobackup/users/ssilvest/perlmutter-test"
# export JULIA_DEPOT_PATH=":${COMMON}/depot"
# export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/satori"
# export JULIA_CUDA_MEMORY_POOL=none
export JULIA=julia

# ====================================================== #
# ================ PREPARATORY ACTIONS ================= #
# ====================================================== #

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
if test -f $BATHYMETRY || test $EXPERIMENT != "RealisticOcean" ; then
   echo "the bathymetry file already exists or we are running the DoubleDrake experiment."
else
    echo "regenerating bathymetry"
    $JULIA --project --check-bounds=no generate_bathymetry.jl
fi

rm generate_bathymetry.jl

# #####
# ##### Now we need to generate daily fluxes and initial conditions
# #####

cd ../data

# check we want to regenerate fluxes otherwise proceed
cat > write_fluxes.jl << EoF_s
include("generate_fluxes.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
generate_fluxes(res; arch = CPU())
EoF_s

if test $REGENERATEFLUXES == 0 ; then
    echo "we are not regenerating the fluxes"
else
    if test $EXPERIMENT != "RealisticOcean" || test $WITHFLUXES == 0 ; then
        echo "WARNING!! we are regenerating the fluxes but they are not needed for this simulation"
    fi
    echo "regenerating fluxes"
    $JULIA --project --check-bounds=no write_fluxes.jl
fi

rm write_fluxes.jl

# check we want to regrid initial conditions otherwise proceed
cat > generate_initial_conditions.jl << EoF_s
include("regrid_initial_conditions.jl")
res = parse(Int, get(ENV, "RESOLUTION", "3"))
Nz  = parse(Int, get(ENV, "NZ", "100"))
regrid_initial_conditions(res, Nz; arch = CPU())
EoF_s

if test $REGENERATEINITIALCONDITIONS == 0 ; then
    echo "we are not regenerating the initial conditions"
else
    if test ! -z $RESTART  || test $EXPERIMENT != "RealisticOcean" ; then
        echo "WARNING!! we are regenerating the initial conditions but it is not needed for this simulation"
    fi
    echo "regenerating initial conditions"
    $JULIA --project --check-bounds=no generate_initial_conditions.jl
fi

rm generate_initial_conditions.jl

# ====================================================== #
# ============ LET'S RUN THE SIMULATION!!! ============= #
# ====================================================== #

#####
##### Now we can finally run our simulation
#####

# cd ../
# sbatch -N ${NNODES} satori_job.sh
