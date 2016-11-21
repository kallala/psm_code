from header import *

nx=64
ny=32*2
nz=32*2
Lx=1.
Ly=1.
Lz=1.
dt=0.000002658*100
c=3.*math.pow(10,8)
c=1
E_0_y=10.
E_0_z=1.
mod=2
w=2*math.pi*c*mod/Lx
vecteur_onde=w/c
print "vecteur d onde %f "%vecteur_onde
S=Simulation(dt,nx,ny,nz,Lx,Ly,Lz,c)
S.Initial_Conditions(1.,E_0_y,E_0_z,mod,0.,0.,0.)

w=2*math.pi*c*mod/S.L_x
K=w/c
phi_y=0.5
phi_z=0.2
K=np.zeros((3,32,32,32))
for i in range(32):
	for j in range(32):
		for k in range(32):
			K[0,j,i,k]=w/c

j=complex(0,1)
#E=fft_champs_3d(S.E_old)
#A=j*ifft_champs_3d(prod_vec_champs_totaux(S.mesh_KK,fft_champs_3d(S.E_old)))
#AA=-j*prod_vec_champs_totaux(K,S.E_old)
nsteps=100
cx=2*math.sin(w*dt/2)
#cx=w*dt
#S.Time_Integration_PSM(nsteps)
for i in range(nsteps):
	if i%100==0:
		print "i= %d " %i
	rotB=j*prod_vec_champs_totaux(S.mesh_KK,S.Btilde_old) #rot B sur le domaine spectral
	#print " verificaiton %f %f "%(np.sum((rotB*S.Btilde_old).real),np.sum((rotB*S.Btilde_old).imag))
	S.Etilde_n=S.Etilde_old+cx*c*c*rotB
	rotE=j*prod_vec_champs_totaux(S.mesh_KK,S.Etilde_n)
	S.Btilde_n=-cx*rotE+S.Btilde_old
	S.Btilde_old=np.copy(S.Btilde_n)
	S.Etilde_old=np.copy(S.Etilde_n)
S.E_n=ifft_champs_3d(S.Etilde_n)
S.B_n=ifft_champs_3d(S.Btilde_n)
x=np.arange(nx)
x=1./nx*x
plt.plot(x,S.E_n[1,:,0,0],"b")
plt.grid(True)
Sa=np.zeros(nx)
for i in range(nx):
	Sa[i]=10.*np.cos(vecteur_onde*x[i]+w*dt*nsteps)
plt.plot(x,Sa,"r")
plt.show()
