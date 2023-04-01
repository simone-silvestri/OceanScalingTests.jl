#!/bin/bash
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH --ntasks-per-node=4
#SBATCH --time=06:00:00
#SBATCH --account=m1759
#SBATCH --gpus-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

module load cray-mpich

export SBATCH_ACCOUNT=m4367
export SALLOC_ACCOUNT=m4367
export COMMON=/global/common/software/m4367

export PATH=${COMMON}/julia-1.9.0-rc2/bin:${PATH}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/perlmutter"
export JULIA_DEPOT_PATH=":${COMMON}/depot"

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

nsys profile --trace=cuda,mpi,nvtx --mpi-impl=mpich \
    srun ./launch.sh julia --check-bounds=no --project experiments/run.jl ${RESOLUTION:=3}
