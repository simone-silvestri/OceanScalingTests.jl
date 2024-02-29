module OceanScalingTests

using Statistics
using JLD2
using Printf
using NVTX
using Oceananigans
using Oceananigans.Units

using Oceananigans.DistributedComputations
using Oceananigans.Fields: interpolate, Field
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.ImmersedBoundaries: PartialCellBottom
using Oceananigans.BoundaryConditions
using CUDA: @allowscalar
using Oceananigans.Operators
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans: prognostic_fields
using MPI

include("correct_oceananigans.jl")
include("bathymetry_and_forcings.jl")
include("boundary_and_initial_conditions.jl")
include("grid_load_balance.jl")
include("quasi_geostrophic_leith.jl")
include("near_global_simulation.jl")
include("simulation_outputs.jl")
include("utils.jl")

end # module OceanScalingTests
