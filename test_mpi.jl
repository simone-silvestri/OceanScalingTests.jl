using CUDA
using MPI
using Statistics: mean

MPI.Init()
using OceanScalingTests

rank = MPI.Comm_rank(MPI.COMM_WORLD)
size = MPI.Comm_size(MPI.COMM_WORLD)

local_comm = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, rank)
node_rank  = MPI.Comm_rank(local_comm)
CUDA.device!(node_rank)

buff = (send = CuArray(zeros(5, 5, 5)), recv = CuArray(zeros(5, 5, 5)))
buff.send .= rank + 10

rank_recv_from = rank - 1 < 0 ? size - 1 : rank - 1
rank_send_to   = rank + 1 > size - 1 ? 0 : rank + 1

recvreq = MPI.Irecv!(buff.recv, rank_recv_from, 0, MPI.COMM_WORLD)
sendreq =  MPI.Isend(buff.send, rank_send_to,   0, MPI.COMM_WORLD)

MPI.Waitall([sendreq, recvreq])

@info rank CUDA.device() mean(buff.send) mean(buff.recv) rank_recv_from rank_send_to
MPI.Barrier(MPI.COMM_WORLD)

MPI.Finalize()
