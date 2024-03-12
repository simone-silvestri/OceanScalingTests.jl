#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500GB
#SBATCH --time 24:00:00

## modules setup
source ../satori/setup_satori.sh

cat > compress.jl << EoF_s
using OceanScalingTests
OceanScalingTests.compress_all_restarts((1080, 450, 100), 8, "./"; prefix = "DoubleDrake_checkpoint_")
OceanScalingTests.compress_surface_fields((1080, 450, 100), 8, "./"; prefix = "DoubleDrake_fields_")
EoF_s

$JULIA --check-bounds=no --project compress.jl
