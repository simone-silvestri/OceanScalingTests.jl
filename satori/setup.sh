module purge all
module add spack
module add cuda/10.1
module load openmpi/3.1.6-cuda-pmi-ucx-slurm-jhklron

export OMPI_MCA_pml=^ucx
export OMPI_MCA_osc=^ucx
export OMPI_MCA_btl_openib_allow_ib=true
