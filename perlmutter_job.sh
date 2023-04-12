#!/bin/bash
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-node=4
#SBATCH --account=m4367
#SBATCH --gpus-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

module load cray-mpich

export SBATCH_ACCOUNT=m4367
export SALLOC_ACCOUNT=m4367

export COMMON=/global/homes/s/ssilvest

export PATH=${COMMON}/julia-1.9-src/bin:${PATH}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/perlmutter"
export JULIA_CUDA_MEMORY_POOL=none

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

##srun --ntasks-per-node=1 dcgmi profile --pause
##srun --ntasks-per-node=1 dcgmi profile --resume

echo "$EXPERIMENT"

export JULIA_GPUCOMPILER_CACHE=$EXPERIMENT
##export JULIA_DEBUG="GPUCompiler"

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

## srun nsys profile --trace=cuda,mpi,nvtx --mpi-impl=mpich --gpu-metrics-device=all \
srun ./launch.sh julia --check-bounds=no --project experiments/run.jl ${RESOLUTION:=3}



