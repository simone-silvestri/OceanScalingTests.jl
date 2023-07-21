using DataDeps

# Initial conditions from 10 years evolved 1/12th degree simulation
path = "https://dl.dropboxusercontent.com/s/4m56hv6hzef5j0j/evolved-initial-conditions-365days.jld2?dl=0"

dh = DataDep("initial_conditions", "1/12th degree resolution initial conditions", path)

DataDeps.register(dh)

datadep"initial_conditions"

# Initial conditions from ECCO version 2 (cube 92) at 01/01/1995
path = "https://dl.dropboxusercontent.com/scl/fi/vq0dz2xo5xnxc0rfzc3fc/ecco-initial-conditions-19950101.jld2?rlkey=kx9ttz728a8unzf40jq2tyzfh&dl=0"
dh = DataDep("ecco_initial_conditions", "1/4th degree resolution initial conditions from ECCO", path)

DataDeps.register(dh)

datadep"ecco_initial_conditions"
