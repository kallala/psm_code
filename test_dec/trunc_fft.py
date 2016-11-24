import numpy as np
import scipy as sp
import math  as math
import cmath as cmath
import matplotlib.pyplot as plt
from mpi4py import MPI
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

nx=128
ngard=16
n_sub=nx/size+ngard
dx=1./nx
if rank==0:
	a=np.arange(nx/size+ngard)
if rank==1:
	a=np.arange(nx/size-ngard,nx,1)
X=np.zeros(n_sub)
for i in range(n_sub):
	X[i]=a[i]*dx
Z=np.arange(nx)

Z=np.cos(12.*dx*np.pi*Z)
Z_loc=Z[np.arange(nx/size*rank,nx/size*(rank+1),1)]
Z_loc=np.cos(12.*math.pi*Z_loc)
z=np.zeros(ngard)
if rank==0:
	Z_loc=np.append(Z_loc,z)
if rank==1:
	Z_loc=np.append(z,Z_loc)


