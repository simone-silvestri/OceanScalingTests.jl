#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time 24:00:00
#SBATCH --reservation=gpu-aware-mpi-testing
#SBATCH --partition=reservation7 
#SBATCH --qos=reservation7 

module purge all
module add spack
module add cuda/11.4
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_NVTX_CALLBACKS=gc
export JULIA_DEPOT_PATH="/nobackup/users/ssilvest/julia_pkg/"
export JULIA_CUDA_MEMORY_POOL=none

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

##nsys profile --trace=nvtx,cuda,mpi --output=report_waitall_$SLURM_JOB_NUM_NODES 

cp satori/* .
srun --mpi=pmi2 ./launch.sh julia --check-bounds=no --project experiments/run.jl ${RESOLUTION:=3}
rm *.toml
