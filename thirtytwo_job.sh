#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH --threads-per-core=1
#SBATCH --mem=1TB
#SBATCH --time 06:00:00
#SBATCH -e error_thirtytwo.txt
#SBATCH -o output_thirtytwo.txt
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

export CUDA_VISIBLE_DEVICES=0,1,2,3

srun -n 64 /nobackup/users/ssilvest/julia-src/julia --check-bounds=no --project experiments/res_thirtytwo.jl
