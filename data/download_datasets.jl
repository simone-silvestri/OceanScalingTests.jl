using DataDeps

# Bathymetry data from ETOPO1 converted to .jld2
path = "https://dl.dropboxusercontent.com/s/4m56hv6hzef5j0j/evolved-initial-conditions-365days.jld2?dl=0"

dh = DataDep("initial_conditions", "1/12th degree resolution initial conditions", path)

DataDeps.register(dh)

datadep"initial_conditions"
