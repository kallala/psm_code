import scipy as sp
import math  as math
import cmath as cmath
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD

def prod_scalaire(k,E,direction):
	nx=E.shape[0]
	a=0*E
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
nx=512
ngard_cells=8
dx=Lx/nx
nsteps=1200
local_coord=np.zeros(nx/(size))
nx_loc_central=nx/size
nx_loc_tot=nx/size+2*ngard_cells

a_loc=max(0,(rank*nx/size-ngard_cells)*dx)
p=(rank*nx/size-ngard_cells)%nx
q=((rank+1)*nx/size+ngard_cells)%nx
local_coord=np.arange(p,q+nx,1)%nx
local_coord=dx*local_coord
#for i in range(local_coord.shape[0]):
#	print "je suis le process %d et coord i = %f"%(rank,local_coord[i])

dt= 2.652582384864922e-05
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
for i in range(ngard_cells,ngard_cells+nx/size,1):
	E_old[i]=cmath.exp(complex(0,-k*(local_coord[i])))
	B_old[i]=cmath.exp(complex(0,-k*(local_coord[i])))

E_old=E0*E_old
#plt.plot(E_old[np.arange(ngard_cells,ngard_cells+nx/size,1)],"g")
#plt.title(rank)
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

Etilde_n=np.copy(Etilde_old)
Btilde_n=np.copy(Btilde_old)

x=range(nx_loc_tot/2)
rx=range(-nx_loc_tot/2+1,0,1)
x.append(0)
x.extend(rx)
local_size=Lx*(1./size+2.*ngard_cells/nx)
K_mesh=2*math.pi/(local_size)*np.asarray(x)

#K_mesh=np.arange(nx)
#K_mesh=np.arange(-nx/2,nx/2,1)
#K_mesh=2*1./Lx*math.pi*K_mesh
j=complex(0,1)
cx=2*math.sin(w*dt/2)
#cx=w*dt
rank_back=(rank-1)%size
rank_front=(rank+1)%size
npa=ngard_cells+nx/size
index_to_send_to_back=np.arange(ngard_cells)
index_to_send_to_front=np.arange(nx/size+ngard_cells,nx/size+2*ngard_cells,1)

for i in range(nsteps):
	if i%1==0:
		print "i= %d " %i
	rotB=j*prod_scalaire(K_mesh,Btilde_n,2)
	Etilde_n=cx/w*(c*c)*rotB+Etilde_n
	rotE=j*prod_scalaire(K_mesh,Etilde_n,1)
	Btilde_n=-cx/w*rotE+Btilde_n
	E_n=np.fft.ifft(Etilde_n)
	B_n=np.fft.ifft(Btilde_n)
	comm.Barrier()
	if i%1==0:
		E_to_send_to_front=E_n[index_to_send_to_front]
		B_to_send_to_front=B_n[index_to_send_to_front]
		Buff_front=np.append(E_to_send_to_front,B_to_send_to_front)
		E_to_send_to_back=E_n[index_to_send_to_back]
		B_to_send_to_back=E_n[index_to_send_to_back]
		Buff_back=np.append(E_to_send_to_back,B_to_send_to_back)
		print"size=%d"%Buff_back.size	
		rcv_f=np.empty(2*ngard_cells,dtype=complex)
                rcv_b=np.empty(2*ngard_cells,dtype=complex)
		print"arrive ici %d"%rank
			
		comm.Send(Buff_front,rank_front,0)
		plt.plot(Buff_front,"r")
		plt.title(rank)
		plt.show()
		print"step %d"%rank

		comm.Barrier()
		print"depasse barriere %d %d"%(rank,rank_back)
	 	comm.Recv(rcv_b,rank_back,0)
		plt.plot(rcv_b)
		plt.title(rank)
		plt.show()
		print"premier rcv %d"%rank	
		comm.Barrier()
		print "depasse 2eme barr %d"%rank
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
			E_n[nx/size+ngard_cells+ind]=0.
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

print"je suis le proc %d"%rank
#Final_loc_E=np.zeros(nx,dtype=complex)
#for i in range(E_n.size):
#	Final_loc_E[nx/size*rank+i]=E_n[i]
#Final_glob_E=np.zeros(nx)
#plt.plot((rank+1)*Final_loc_E)
#plt.show()
#comm.Reduce(Final_loc_E,Final_glob_E,op=MPI.SUM,root=0)
if rank==0:
	E_field=np.zeros(nx,dtype=complex)
comm.Barrier()

if rank==0:
	for k in range(nx/size):
			for h in range(size):
				E_field[k+h*nx/size]=data2[h][k]
	#	E_field[k]=E_n[k]
	#	E_field[k+2*nx/size]=A[k]
	#	E_field[k+nx/size]=B[k]
	#	E_field[k+3*nx/size]=C[k]
#Diag=np.zeros((nx,nx))  
#Solution=np.zeros(nx)
#phase=w*nsteps*dt
k=w/c
if rank == 0:
	Solution=np.zeros(nx,dtype=complex)
	phase=w*dt*nsteps
	for i in range (nx):
		Solution[i]=E0*math.cos(-phase-k*X[i])
#for i in range(nx):     
#	Diag[i,i]=K_mesh[i]
#Diag=np.asmatrix(Diag)

#print "phase initiale exacte %f"%math.acos(math.cos(phase))
#print "phase initiale simulee %f"%math.acos((1./E0*E_n[0]).real)
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












