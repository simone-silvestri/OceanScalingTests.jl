__precompile__()
module GenerateBathymetry

export interpolate_bathymetry_from_ETOPO1, write_bathymetry_to_file
export LinearInterpolation, SplineInterpolation, SpectralInterpolation

using DataDeps
using PyCall
using JLD2
using Oceananigans
using Oceananigans.Fields: interpolate, set!
using Oceananigans.Architectures: device, arch_array
using Oceananigans.Units
using Oceananigans.Fields: interpolate
using Statistics: dot
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.Grids: architecture

# a small number larger than 0, 
# to avoid finite precision errors
const ABOVE_SEA_LEVEL = 0.001

include("utils.jl")
include("interpolate_bathymetry.jl")

const sckikitimage = PyNULL()

function __init__(; remove_existing_data=false)
    copy!(sckikitimage, pyimport_conda("skimage.measure", "scikit-image"))

    # Bathymetry data from ETOPO1 converted to .jld2
    path = "https://dl.dropboxusercontent.com/s/33v6tye7sftt1ud/bathymetry-ice-21600x10800.jld2?dl=0"

    dh = DataDep("etopo1_bathymetry", "1/60th degree resolution bathymetry", path)

    DataDeps.register(dh)

    datadep"etopo1_bathymetry"

    remove_existing_data && rm(datadep"etopo1_bathymetry", recursive=true, force=true)
end

end # module
