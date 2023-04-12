#!/bin/bash

export RESOLUTION=48
export LOADBALANCE=1
export NZ=100
export EXPERIMENT="RealisticOcean"
export PROFILE=1
export WITHFLUXES=0

sbatch -N144 perlmutter_job.sh
