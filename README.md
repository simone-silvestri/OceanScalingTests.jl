1. Add Julia to path export PATH=/TO/JULIA:${PATH}
2. EXPERIMENT: 
    - "Quiescent" no flow, no bathymetry and zero flux boundary conditions
    - "DoubleDrake" setup similar to [https://doi.org/10.1175/2009JCLI3197.1](https://doi.org/10.1175/2009JCLI3197.1), default
    - "RealisticOcean" a realistic ocean setup with initial conditions and boundary forcing from ECCO2 Version 4 climatological data. It assumes data is available to load in the `data` folder
3. PROFILE: 0 (run simulation), 1 (profile 100 time steps marked with NVTX, default) 
4. NZ: vertical levels (default 120)
5. WITHFLUXES: 0, 1 (only relevant for RealisticOcean, default 1)
6. RESTART: restart file iterations (only relevant for RealisticOcean, default "")

## Satori
`EXPERIMENT="DoubleDrake" RESOLUTION=6 sbatch -N1 satori_job.sh`

## Packages:
- `Oceananigans.jl.git#ss-vc/overlay`
- `CUDA.jl.git#vc/ka_transition`

## Example

vorticity (left) and free surface (right) for the "DoubleDrake" experiment at resolution = 3 ran on 2 mpi processors (9 years integration)


https://user-images.githubusercontent.com/33547697/219111104-666e1215-933d-4d31-9628-57ed87b335e4.mp4

