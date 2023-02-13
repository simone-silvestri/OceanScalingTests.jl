#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --threads-per-core=1
#SBATCH --mem=200GB
#SBATCH --time 06:00:00
#SBATCH --account=m1759

module load cuda
module load cray-mpich
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

export SBATCH_ACCOUNT=m1759
export SALLOC_ACCOUNT=m1759
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

cp perlmutter/* .
srun --mpi=pmi2 --gpus ./launch.sh julia --check-bounds=no --project= experiments/run.jl ${RESOLUTION:=3}
rm *.toml
