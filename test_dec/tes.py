from mpi4py import MPI
import numpy as np
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
data =np.zeros(5)
data =data+ (rank+1)**2
data2 = comm.gather(data, root=0)
if rank == 0:
	print"size to %d"%data2[0][0]
	for i in range(size):
		for j in range(5):
			print"hel %d"%data2[i][j]
			#print" i= %d data= %d"%(i, data2[i,j])
