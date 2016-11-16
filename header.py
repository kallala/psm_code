from tools import *

class Simulation :  # class simulation solveur EM 1D pseudo-specral



	def   __init__(self,dt,nx,ny,nz,Lx,Ly,Lz):
		c =  3*math.pow(10,8)
		self.n_x     = nx
		self.n_y     = ny
		self.n_z     = nz
		self.L_x     = Lx
		self.L_y     = Ly
		self.L_z     = Lz
		self.d_t     = dt
 
		self.E_n     = np.zeros((nx,ny,nz,3),dtype=complex)
		self.Etilde_n= np.zeros((nx,ny,nz,3),dtype=complex)

		self.E_old   = np.zeros((nx,ny,nz,3),dtype=complex)
		self.Etilde_old = np.zeros((nx,ny,nz,3),dtype=complex)

		self.B_n     = np.zeros((nx,ny,nz,3),dtype=complex)
		self.Btilde_n= np.zeros((nx,ny,nz,3),dtype=complex)
	
		self.B_old   = np.zeros((nx,ny,nz,3),dtype=complex)
		self.Btilde_old = np.zeros((nx,ny,nz,3),dtype=complex)	
		self.d_x     = (Lx-0)/nx
	   	self.d_y     = (Ly-0)/ny
		self.d_z     = (Lz-0)/nz
		#verification de la cfl	
		alpha=2/math.pi*math.sqrt(1./(self.d_x*self.d_x)+1./(self.d_y*self.d_y)+1./(self.d_z*self.d_z))/self.d_t
		if alpha>c: 
			print " la condition cfl est verifiee cfl  %f " %(c/alpha)
		else:
			print " la condition cfl est non verifiee cfl = %f " %(c/alpha)
		#construct physical mesh
		
		X=range(self.n_x)
		X=np.asarray(X)
		X=self.d_x*X

		Y=range(self.n_y)
		Y=np.asarray(Y)
		Y=self.d_y*Y

		Z=range(self.n_z)
		Z=np.asarray(Z)
		Z=self.d_z*Z

		self.mesh=np.meshgrid(X,Y,Z)
	
		
		#construct spectral mesh
		X=range(self.n_x/2)
		rX=range(-self.n_x/2+1,0,1)
		X.append(0)
		X.extend(rX)
		X=np.asarray(X)
		X=2*math.pi/(self.L_x-0)*X
		Y=range(self.n_y/2)	
		rY=range(-self.n_y/2+1,0,1)
		Y.append(0)
		Y.extend(rY)
		Y=np.asarray(Y)
		Y=2*math.pi/(self.L_y-0)*Y
		Z=range(self.n_z/2)	
		rZ=range(-self.n_z/2+1,0,1)
		Z.append(0)
		Z.extend(rZ)
		Z=np.asarray(Z)
		Z=2*math.pi/(self.L_z-0)*Z
	        self.mesh_k=np.meshgrid(X,Y,Z)
	def Initial_Conditions(self,E_0_x,E_0_y,E_0_z,mod,phi_x,phi_y,phi_z):
		c=3.*math.pow(10,8)
		w=2*math.pi*c*mod/self.L_x
		vecteur_onde=w/c
		for i in range(self.n_x):
			for j in range(self.n_y):
				for k in range( self.n_z):
					E1=E_0_y*cmath.exp(complex(0,-vecteur_onde*self.mesh[0][j,i,k]+phi_y))
					self.E_old[i,j,k,1]=E1
					E2=E_0_z*cmath.exp(complex(0,-vecteur_onde*self.mesh[0][i,i,k]+phi_z))
					self.E_old[i,j,k,2]=E2
		Etilde_old=fft_champs_3d(self.E_old)
		Btilde_old=fft_champs_3d(self.B_old)
		

	def Time_Integration_PSM(self):
		print "test %d " %self.n_x		
