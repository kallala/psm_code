from mpi4py import MPI
import numpy as np
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size=comm.Get_size()
data=np.zeros(4)+rank

dest=(rank+1)%size
comm.send(data, dest, tag=11)
print"le proc %d a envoy"%rank
comm.Barrier()
source=(rank-1)%size
data = comm.recv(source, tag=11)
for i in range(data.size):
	print("recu %d mon rang %d ") %(data[i],rank)
