# using Preferences
# const iscray = parse(Bool, load_preference(Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"), "iscray", "false"))
# @debug "Preloading GTL library" iscray
# if iscray
#     import Libdl
#     Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
# end
# using OceanScalingTests
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using NVTX
using MPI
using JLD2
using CUDA 

using Oceananigans.Models.HydrostaticFreeSurfaceModels

MPI.Init()

comm   = MPI.COMM_WORLD
rank   = MPI.Comm_rank(comm)
Nranks = MPI.Comm_size(comm)

Ry = 1
Rx = Nranks

@show Rx, Ry

ranks       = (Rx, Ry, 1)

# Enviromental variables
resolution  = parse(Int,  get(ENV, "RESOLUTION", "3"))
experiment  = Symbol(     get(ENV, "EXPERIMENT", "DoubleDrake"))
with_fluxes = parse(Bool, get(ENV, "WITHFLUXES", "0"))
profile     = parse(Bool, get(ENV, "PROFILE", "1"))
restart     =             get(ENV, "RESTART", "")
Nz          = parse(Int,  get(ENV, "NZ", "100"))
loadbalance = parse(Bool, get(ENV, "LOADBALANCE", "0"))
precision   = eval(Symbol(get(ENV, "PRECISION", "Float64")))

using Oceananigans.Architectures: arch_array
using Oceananigans.Distributed

buff = (send = arch_array(CPU(), zeros(5, 100, 110)), recv = arch_array(CPU(), zeros(5, 100, 110)))

arch = DistributedArch(CPU(); ranks, topology = (Periodic, Bounded, Bounded))

rank_to_send_to   = arch.connectivity.west
rank_to_recv_from = arch.connectivity.east

@info MPI.Comm_rank(MPI.COMM_WORLD) CUDA.device() rank_to_send_to rank_to_recv_from

buff.send .= MPI.Comm_rank(MPI.COMM_WORLD)

reqrecv = MPI.Irecv!(buff.recv, rank_to_recv_from, 0, arch.communicator)
reqsend =  MPI.Isend(buff.send, rank_to_send_to,   0, arch.communicator)

MPI.Waitall([reqrecv, reqsend])

using Statistics: mean

@info MPI.Comm_rank(MPI.COMM_WORLD) mean(buff.send) mean(buff.recv)

MPI.Barrier(MPI.COMM_WORLD)

MPI.Finalize()
