#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --time 24:00:00
#SBATCH --reservation=gpu-aware-mpi-testing
#SBATCH --partition=reservation7 
#SBATCH --qos=reservation7 
###SBATCH -w node[0001-0003,0006,0010,0013-0019,0021-0028,0030-0033,0036-0037]

module purge all
module add spack
module add cuda/11.4
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

export COMMON="/nobackup/users/ssilvest/perlmutter-test"

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_NVTX_CALLBACKS=gc
export JULIA_DEPOT_PATH=":${COMMON}/depot"
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/satori"
export JULIA_CUDA_MEMORY_POOL=none
export JULIA="${COMMON}/julia-1.9-src/julia"

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

## 

nsys profile --trace=nvtx,cuda,mpi --output=${COMMON}/report_N${SLURM_JOB_NUM_NODES}_R${RESOLUTION}_${PRECISION} srun --mpi=pmi2 ./launch.sh $JULIA --check-bounds=no --project experiments/run.jl ${RESOLUTION:=3}
