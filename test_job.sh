#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --mem=100GB
#SBATCH --time 24:00:00

## modules setup
source satori/setup_satori.sh

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

$NSYS srun --mpi=pmi2 ./launch.sh $JULIA --check-bounds=no --project test_mpi.jl ${RESOLUTION:=3}
