1. Add Julia to path export PATH=/TO/JULIA:${PATH}
2. EXPERIMENT: 
    - "Quiescent" no flow, no bathymetry and zero flux boundary conditions, default
    - "DoubleDrake" setup similar to [https://doi.org/10.1175/2009JCLI3197.1](https://doi.org/10.1175/2009JCLI3197.1)
3. USEBUFFERS: 0 (sending `views`), 1 (buffered communication, default)

## Satori
4. `EXPERIMENT="DoubleDrake" USEBUFFERS=1 RESOLUTION=6 sbatch -N1 satori_job.sh`
## Perlmutter
4. `EXPERIMENT="DoubleDrake" USEBUFFERS=1 RESOLUTION=6 sbatch -N1 perlmutter_job.sh`
