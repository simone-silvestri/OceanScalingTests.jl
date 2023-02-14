#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --threads-per-core=1
#SBATCH --mem=1TB
#SBATCH --time 06:00:00
#SBATCH --reservation=gpu-aware-mpi-testing
#SBATCH --partition=reservation7 
#SBATCH --qos=reservation7 

module purge all
module add spack
module add cuda/10.1
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_NVTX_CALLBACKS=gc

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh


cp satori/* .
nsys --profile --trace=nvtx,cuda,mpi --output=report_new srun --mpi=pmi2 ./launch.sh julia --check-bounds=no --project= experiments/run.jl ${RESOLUTION:=3}
rm *.toml
