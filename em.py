import scipy
import math
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD


def prod_vect(k,E,direction):
	nx=E.shape[0]
	a=0*E
	if(direction==1):#ex^ey
		for i in range(nx):
			a[i]=k[i]*E[i]
	if (direction==2):#ex^ez
		for i in range(nx):
			a[i]=-k[i]*E[i]
	return(a)

# User defined parameters ##################

c  = 1.                                           #Speed of light
Lx = 1.                                           #Domain size
nx = 512                                          #Number of cells = number of points because of periodic boundaries
ngard_cells = 32                                  #Number of guard cell
dx = Lx/nx                                        #Cell size
dt = 2.652582384864922e-05                        #Time step
nsteps = 1200                                     #Number of iterations
nx_loc_tot = nx/size + 2*ngard_cells              #Number of cells per MPI domain

########################################################################################################

print "This program works only if the number of proc ", size, " divides the number of cells ", size

local_coord = scipy.arange(nx_loc_tot)        #Coordinates of the nx_loc_tot local MPI points.
local_coord += rank*nx/size - ngard_cells     #nx/size non guard points per MPI and the first non guard cell is located at x=0.
local_coord = (local_coord % nx) * dx         #periodic boundaries and account for individual cell size dx. 


if (rank==0):
	cfl=c*dt/dx
	print"cfl= %f "%cfl
X=np.arange(nx)*Lx*dx
mode=3
w=2*math.pi*c*mode/Lx
k=w/c
E0=100.
B0=-E0/c
E_old=np.zeros(nx_loc_tot,dtype=complex)
B_old=0.*E_old



for i in range(ngard_cells,ngard_cells+nx/size,1):
	E_old[i]=scipy.exp(-1j*k*(local_coord[i]))
	B_old[i]=scipy.exp(-1j*k*(local_coord[i]))

E_old=E0*E_old

#Plot initial conditions
#plt.plot(E_old)
#plt.title(rank)
#plt.show()
#plt.plot(E_old[np.arange(ngard_cells,ngard_cells+nx/size,1)],"g")
#plt.title(rank)

B_old=B0*B_old
E_n=0*E_old
B_n=0*B_old

Etilde_old=np.fft.fft(E_old)
Btilde_old=np.fft.fft(B_old)

Etilde_n=np.copy(Etilde_old)
Btilde_n=np.copy(Btilde_old)

x=range(nx_loc_tot/2)
rx=range(-nx_loc_tot/2+1,0,1)
x.append(0)
x.extend(rx)
local_size=Lx*(1./size+2.*ngard_cells/nx)
K_mesh=2*math.pi/(local_size)*np.asarray(x)

cx=2*math.sin(w*dt/2) #PSATD
#cx=w*dt #PSTD

rank_back=(rank-1)%size
rank_front=(rank+1)%size
npa=ngard_cells+nx/size
index_to_send_to_back=np.arange(ngard_cells)
index_to_send_to_front=np.arange(ngard_cells)+ngard_cells+nx/size

for i in range(nsteps):
	if (i%100==0 and rank ==0) :
		print "iteration = %d " %i
	rotB=1j*prod_vect(K_mesh,Btilde_n,2)
	Etilde_n=cx/w*(c*c)*rotB+Etilde_n
	rotE=1j*prod_vect(K_mesh,Etilde_n,1)
	Btilde_n=-cx/w*rotE+Btilde_n
	E_n=np.fft.ifft(Etilde_n)
	B_n=np.fft.ifft(Btilde_n)
	comm.Barrier()
	if i%1==0:
		E_to_send_to_front=E_n[index_to_send_to_front]
		B_to_send_to_front=B_n[index_to_send_to_front]
		Buff_front=np.append(E_to_send_to_front,B_to_send_to_front)
		E_to_send_to_back=E_n[index_to_send_to_back]
		B_to_send_to_back=B_n[index_to_send_to_back]
		Buff_back=np.append(E_to_send_to_back,B_to_send_to_back)
		rcv_f=np.empty(2*ngard_cells,dtype=complex)
                rcv_b=np.empty(2*ngard_cells,dtype=complex)
		comm.Send(Buff_front,rank_front,0)
		comm.Barrier()
	 	comm.Recv(rcv_b,rank_back,0)
		comm.Barrier()
		comm.Send(Buff_back,rank_back,0)
		comm.Barrier()
		comm.Recv(rcv_f,rank_front,0)
		comm.Barrier()
		for ind in range(ngard_cells):
			E_n[ngard_cells+ind]+=rcv_b[ind]
			B_n[ngard_cells+ind]+=rcv_b[ngard_cells+ind]
			E_n[nx/size+ind]+=rcv_f[ind]
			B_n[nx/size+ind]+=rcv_f[ngard_cells+ind]
			E_n[ind]=0.
			B_n[ind]=0.
			E_n[nx/size+ngard_cells+ind]=0.
			B_n[nx/size+ngard_cells+ind]=0.
		Etilde_n=np.fft.fft(E_n)
		Btilde_n=np.fft.fft(B_n)
	comm.Barrier()
true_coord=np.arange(nx/size)+ngard_cells
E_n=E_n[true_coord]
plt.plot(E_n)
plt.title(rank)
plt.grid(True)
plt.show()
B_n=B_n[true_coord]

data2=comm.gather(E_n,0)

if rank==0:
	E_field=np.zeros(nx,dtype=complex)
comm.Barrier()

if rank==0:
	for k in range(nx/size):
			for h in range(size):
				E_field[k+h*nx/size]=data2[h][k]
k=w/c
if rank == 0:
	Solution=np.zeros(nx,dtype=complex)
	phase=w*dt*nsteps
	for i in range (nx):
		Solution[i]=E0*math.cos(-phase-k*X[i])
	plt.plot(Solution,'bo')
comm.Barrier()
if rank ==0:
	plt.grid(True)
	plt.plot(E_field,'r')
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












