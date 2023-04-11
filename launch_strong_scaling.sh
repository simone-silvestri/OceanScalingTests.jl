export EXPERIMENT="DoubleDrake"
export PROFILE=1
export WITHFLUXES=0

export LOADBALANCE=0

FOLDER="reports" RESOLUTION=12 PRECISION="Float32" sbatch -N1 satori_cached_job.sh
