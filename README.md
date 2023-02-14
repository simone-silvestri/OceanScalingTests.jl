1. Add Julia to path export PATH=/TO/JULIA:${PATH}
2. EXPERIMENT: 
    - "Quiescent" no flow, no bathymetry and zero flux boundary conditions, default
    - "DoubleDrake" setup similar to [https://doi.org/10.1175/2009JCLI3197.1](https://doi.org/10.1175/2009JCLI3197.1)
    - "RealisticOcean" a realistic ocean setup with initial conditions and boundary forcing from ECCO2 Version 4 climatological data. It assumes data is available to load in the `data` folder
3. USEBUFFERS: 0 (sending `views`), 1 (buffered communication, default)

## Satori
4. `EXPERIMENT="DoubleDrake" RESOLUTION=6 sbatch -N1 satori_job.sh`
## Perlmutter
4. `EXPERIMENT="DoubleDrake" RESOLUTION=6 sbatch -N1 perlmutter_job.sh`

## Example

vorticity (left) and free surface (right) for the "DoubleDrake" experiment at resolution = 3 ran on 2 mpi processors (5 years integration)


https://user-images.githubusercontent.com/33547697/218850549-f4efc360-7937-4833-8aef-c4dd4863a605.mp4

