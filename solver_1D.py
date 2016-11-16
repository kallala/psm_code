from header import *

nx=128
ny=128
nz=128
Lx=1.
Ly=1.
Lz=1.
dt=0.00001
c=3.*math.pow(10,8)
E_0_y=1.
E_0_z=1.
S=Simulation(dt,nx,ny,nz,Lx,Ly,Lz)
mod = 4
w=2*math.pi*c*mod/S.L_x
K=w/c
phi_y=0.5
phi_z=0.2



for i in range(nx):
	for j in range (ny):
		for k in range (nz):
			E1=E_0_y*cmath.exp(complex(0,-k*S.mesh_x[i,j,k])+phi_y)
			S.E_old[i,j,k,1]=E1
			E2=E_0_z*cmath.exp(complex(0,-k*S.mesh_x[i,j,k])+phi_z)
			S.E_old[i,j,k,2]=E2
ml=S.Time_Integration_PSM()
			 





