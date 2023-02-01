module OceanScalingTests

using Statistics
using JLD2
using Printf
using Oceananigans
using Oceananigans.Units

using Oceananigans.Distributed
using Oceananigans.Fields: interpolate, Field
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.BoundaryConditions
using CUDA: @allowscalar
using Oceananigans.Operators
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans: prognostic_fields
using MPI

include("vertical_faces.jl")
include("near_global_simulation.jl")

end # module OceanScalingTests
