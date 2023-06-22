module load cray-mpich

export COMMON=/global/homes/s/ssilvest

export PATH=${COMMON}/julia-1.9-src/bin:${PATH}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/perlmutter"

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

export RESOLUTION=1
export NZ=12
export PROFILE=1

echo "$EXPERIMENT"

export JULIA_GPUCOMPILER_CACHE=$EXPERIMENT

julia --check-bounds=no --project -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"

rm -rf ${COMMON}/.julia/scratchspaces

julia --check-bounds=no --project experiments/run.jl
