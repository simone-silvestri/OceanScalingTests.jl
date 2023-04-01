module load cray-mpich

export COMMON=/global/common/software/m4367

export PATH=${COMMON}/julia-1.9.0-rc2/bin:${PATH}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_LOAD_PATH="${JULIA_LOAD_PATH}:$(pwd)/perlmutter"
