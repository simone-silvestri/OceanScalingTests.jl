# Satori

1. Add Julia to path export PATH=/TO/JULIA:${PATH}
2. EXPERIMENT: 
    - "Quiescent" (no flow, no bathymetry, zero flux boundary conditions)
    - "DoubleDrake" setup similar to [https://doi.org/10.1175/2009JCLI3197.1](https://doi.org/10.1175/2009JCLI3197.1)
3. USEBUFFERS: 0 (sending `views`), 1 (buffered communication)
4. `EXPERIMENT="DoubleDrake" USEBUFFERS=0 RESOLUTION=6 sbatch -N1 satori_job.sh`