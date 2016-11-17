from header import *

nx=128
ny=128
nz=128
Lx=1.
Ly=1.
Lz=1.
dt=00.000000015
c=3.*math.pow(10,8)
E_0_y=1.
E_0_z=1.
mod=2
w=2*math.pi*c*mod/Lx
vecteur_onde=w/c
S=Simulation(dt,nx,ny,nz,Lx,Ly,Lz)
S.Initial_Conditions(1.,E_0_y,E_0_z,mod,0.,0.,0.)

w=2*math.pi*c*mod/S.L_x
K=w/c
phi_y=0.5
phi_z=0.2



#E=fft_champs_3d(S.E_old)



