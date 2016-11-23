import numpy as np
import scipy as sp
import math  as math
import cmath as cmath
import matplotlib.pyplot as plt
from mpi4py import MPI
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()


def prod_scalaire(k,E,direction,nx):
	a=np.zeros(nx,dtype=complex)
	if(direction==1):#ex^ey
		for i in range(nx):
			a[i]=k[i]*E[i]
	if (direction==2):#ex^ez
		for i in range(nx):
			a[i]=-k[i]*E[i]
#	a=np.asmatrix(a)
#	a=a.T
	return(a)
c=3*math.pow(10,8)
c=1.
Lx=1.
nx=1024
ngard_cells=32
dx=Lx/nx
nsteps=150
local_coord=np.zeros(nx/(size))
nx_loc_central=nx/size
nx_loc_tot=nx/size+2*ngard_cells
if rank ==0:
	nx_loc_tot-=ngard_cells
if rank ==size-1:
	nx_loc_tot-=ngard_cells
a_loc=max(0,(rank*nx/size-ngard_cells)*dx)
local_index=range(max(rank*nx/size-ngard_cells,0),min((rank+1)/size+ngard_cells,nx),1)
for i in local_index:
	local_coord[i]=a_loc+dx*i
dt= 2.652582384864922e-05
if (rank==0):
	cfl=c*dt/dx
	print"cfl= %f "%cfl
X=np.arange(nx)*Lx*dx
mode=7
w=2*math.pi*c*mode/Lx
k=w/c
E0=100.
B0=-E0/c
E_old=np.zeros(nx_loc_tot,dtype=complex)
B_old=0*E_old
#Tf and Tf-1 matrix
#Tf=np.zeros((nx,nx),dtype=complex)
#Tf1=np.zeros((nx,nx),dtype=complex)
#for i in range(nx):
#	for j in range(nx):
#		Tf[i,j]=cmath.exp(complex(0,-2*math.pi*i*j/float(nx)))
#		Tf1[i,j]=1./nx*cmath.exp(complex(0,2*math.pi*i*j/float(nx)))
#Tf=np.asmatrix(Tf)
#Tf1=np.asmatrix(Tf1)
for i in local_index:
	E_old[i]=cmath.exp(complex(0,-k*(a_loc+i*dx)))
	B_old[i]=cmath.exp(complex(0,-k*(a_loc+i*dx)))

E_old=E0*E_old

B_old=B0*B_old
E_n=0*E_old
B_n=0*B_old
#Etilde_old=np.fft.fft(E_old)
#Btilde_old=np.fft.fft(B_old)

#E_old=np.asmatrix(E_old).T
#B_old=np.asmatrix(B_old).T
#E_n=np.asmatrix(E_n).T
#B_n=np.asmatrix(B_n).T

Etilde_old=np.fft.fft(E_old)
Btilde_old=np.fft.fft(B_old)

Etilde_n=0*Etilde_old
Btilde_n=0*Btilde_old

x=range(nx_loc_tot/2)
rx=range(-nx_loc_tot/2+1,0,1)
x.append(0)
x.extend(rx)
local_size=Lx/size+2*ngard_cells*dx
if rank==0:
	local_size=local_size-ngard_cells*dx
if rank==size-1:
	local_size=local_size-ngard_cells*dx
	
K_mesh=2*math.pi/(local_size)*np.asarray(x)

#K_mesh=np.arange(nx)
#K_mesh=np.arange(-nx/2,nx/2,1)
#K_mesh=2*1./Lx*math.pi*K_mesh
j=complex(0,1)
cx=2*math.sin(w*dt/2)
cx=w*dt
for i in range(nsteps):
	if i%1000==0:
		print "i= %d " %i
	rotB=j*prod_scalaire(K_mesh,Btilde_old,2,nx_loc_tot)
	Etilde_n=cx/w*(c*c)*rotB+Etilde_old
	rotE=j*prod_scalaire(K_mesh,Etilde_n,1,nx_loc_tot)
	Btilde_n=-cx/w*rotE+Btilde_old
	Etilde_old=Etilde_n
	Btilde_old=Btilde_n
E_n=np.fft.ifft(Etilde_n)
B_n=np.fft.ifft(Btilde_n)

#Diag=np.zeros((nx,nx))  
Solution=np.zeros(nx)
phase=w*nsteps*dt
for i in range (nx):
	Solution[i]=E0*math.cos(-phase-k*X[i])
#for i in range(nx):     
#	Diag[i,i]=K_mesh[i]
#Diag=np.asmatrix(Diag)

print "phase initiale exacte %f"%math.acos(math.cos(phase))
print "phase initiale simulee %f"%math.acos((1./E0*E_n[0]).real)
plt.grid(True)
#plt.plot(Solution,'bo')
plt.plot(E_n,'r')
plt.show()
#resolution des eq de maxwell en 1D 
#l onde se propage suivant x 
#E est porte par y 
#B est porte par Z
# la seule variable spatiale importante est x 
#rot E=-dB/dt
#rotB=1/(c*c)*dE/dt
#les conditions initiales doivent verifier divE=0 divB=0
#le champs initial est donne par E0exp(-kx)*ey
#B0 donne par B0exp(-kx)ez
#k=w/c
#on a B0=-E/c
#le domaine est de longueur 1












