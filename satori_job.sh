#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=16
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

if $PROFILE; then
   NSYS="nsys profile --trace=nvtx,cuda,mpi --output=${COMMON}/report_N${SLURM_JOB_NUM_NODES}_R${RESOLUTION}_${PRECISION}"
fi

$NSYS srun --mpi=pmi2 ./launch.sh $JULIA --check-bounds=no --project experiments/run.jl ${RESOLUTION:=3}
