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

include("bathymetry_and_forcings.jl")
include("boundary_and_initial_conditions.jl")
include("grid_load_balance.jl")
include("quasi_geostrophic_leith.jl")
include("near_global_simulation.jl")
include("simulation_outputs.jl")
include("utils.jl")

using PrecompileTools
@setup_workload let
    ranks       = (1, 1, 1)

    # Enviromental variables
    resolution     = 1
    experiment     = :DoubleDrake
    with_fluxes    = false
    with_restoring = false
    profile        = false
    restart        = ""
    Nz             = 1
    loadbalance    = false
    precision      = Float64
    output_dir     = mktempdir()

    final_year  = 0
    final_month = 1

    max_Δt = precision(45 * 48 * 1.5 / resolution)
    min_Δt = 5 # precision(45 * 48 * 0.9 / resolution)
    stop_time = 7300days
    @compile_workload begin
        simulation = OceanScalingTests.scaling_test_simulation(resolution, ranks, (min_Δt, max_Δt), stop_time; Nz, experiment, restart,
                                    profile, with_fluxes, with_restoring, loadbalance, precision)

        OceanScalingTests.set_outputs!(simulation, Val(experiment); output_dir, overwrite_existing = false, checkpoint_time = 10days)
        run!(simulation, pickup = false)
    end
end


end # module OceanScalingTests
