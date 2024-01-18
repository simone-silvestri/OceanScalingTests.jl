#!/bin/bash

# Let's first specify the enviromental variables

# ====================================================== #
# ================ USER SPECIFIED INPUTS =============== #
# ====================================================== #

# Grid size
# (RESOLUTION * 360 in x and RESOLUTION * 150 in y and Nz in z)
export RESOLUTION=3 
export NZ=50

# Experiment details
export EXPERIMENT="DoubleDrake"
export PRECISION="Float64"

# Might increase performance when using a lot of cores (i.e. improves scalability)
# Matters only for the RealisticOcean experiment
export LOADBALANCE=0

# Do we want to profile or run a simulation?
export PROFILE=1

# How many nodes are we running on?
export NNODES=2

# Restart from interpolated fields ? "" : "numer_of_iteration_to_restart_from"
export RESTART=""

# Server specific enviromental variables and modules
# Edit depending on the system 
# source satori/setup_satori.sh
export JULIA_CUDA_MEMORY_POOL=none
export JULIA=julia

# Profile specific variable
export JULIA_NVTX_CALLBACKS=gc

# ====================================================== #
# ============ LET'S RUN THE SIMULATION!!! ============= #
# ====================================================== #

#####
##### Now we can finally run our simulation
#####

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=3
exec \$*
EoF_s
chmod +x launch.sh

NSYS="/storage6/simone/new_nsight/bin/nsys"
/storage6/simone/new_nsight/bin/nsys profile --trace=nvtx,cuda --output=./report_test ./launch.sh julia --project --check-bounds=no experiments/run.jl
