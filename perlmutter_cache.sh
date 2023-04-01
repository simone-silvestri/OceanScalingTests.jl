module load cray-mpich

export COMMON=/global/common/software/m4367

export PATH=${COMMON}/julia-1.9.0-rc2/bin:${PATH}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/perlmutter"

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

export RESOLUTION=1
export NZ=12
export EXPERIMENT=DoubleDrake

JULIA_DEPOT_PATH="${COMMON}/depot" julia --check-bounds=no --project -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"

rm -rf ${HOME}/.julia/{packages, compiled}

export JULIA_DEPOT_PATH=":${COMMON}/depot"
julia --check-bounds=no --project experiments/run.jl