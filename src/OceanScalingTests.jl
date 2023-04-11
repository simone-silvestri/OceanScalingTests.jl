module OceanScalingTests

using Statistics
using JLD2
using Printf
using NVTX
using Oceananigans
using Oceananigans.Units

using Oceananigans.Distributed
using Oceananigans.Fields: interpolate, Field
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using Oceananigans.BoundaryConditions
using CUDA: @allowscalar
using Oceananigans.Operators
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans: prognostic_fields
using MPI

include("bathymetry_and_forcings.jl")
include("boundary_and_initial_conditions.jl")
include("grid_load_balance.jl")
include("near_global_simulation.jl")
include("simulation_outputs.jl")

end # module OceanScalingTests
