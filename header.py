from tools import *

class Simulation :  # class simulation solveur EM 1D pseudo-specral



	def   __init__(self,dt,nx,ny,nz,Lx,Ly,Lz):

		self.n_x     = nx
		self.n_y     = ny
		self.n_z     = nz
		self.L_x     = Lx
		self.L_y     = Ly
		self.L_z     = Lz
		self.d_t     = dt
 
		self.E_n     = np.zeros((nx,ny,nz,3),dtype=complex)
	
		self.E_old   = np.zeros((nx,ny,nz,3),dtype=complex)

		self.B_n     = np.zeros((nx,ny,nz,3),dtype=complex)

		self.B_old   = np.zeros((nx,ny,nz,3),dtype=complex)

	        self.mesh    = np.zeros((nx,ny,nz,3))
		
		self.mesh_k  = np.zeros((nx,ny,nz,3))
		
		self.d_x     = (Lx-0)/nx
	   	self.d_y     = (Ly-0)/ny
		self.d_z     = (Lz-0)/nz

	def Construct_Mesh(self):	
		#construct physical mesh
		for i in range (self.n_x):
			for j in range (self.n_y):
				for k in range (nz):
					self.mesh[i,j,k,0] = i*self.d_x
					self.mesh[i,j,k,1] = j*self.d_y
					self.mesh[i,j,k,2] = k*self.d_z
		#construct spectral mesh
		X=range(self.n_x/2)
		rX=range(-self.n_x/2+1,0,1)
		X.append(0)
		X.extend(rX)
		X=np.asarray(X)
		X=2*math.pi/(self.L_x-0)*X
		self.k[:,:,:,0]=X
		X=range(self.n_y/2)	
		rX=range(-self.n_y/2+1,0,1)
		X.append(0)
		X.extend(rX)
		X=np.asarray(X)
		X=2*math.pi/(self.L_y-0)*X
		self.k[:,:,:,1]=X
		X=range(self.n_z/2)	
		rX=range(-self.n_z/2+1,0,1)
		X.append(0)
		X.extend(rX)
		X=np.asarray(X)
		X=2*math.pi/(self.L_z-0)*X
		self.k[:,:,:,2]=X
		





	def Time_Integration_PSM(self):
		print "test %d " %self.n_x		

