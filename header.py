from tools import *

class Simulation :  # class simulation solveur EM 1D pseudo-specral



	def   __init__(self,dt,nx,ny,nz,Lx,Ly,Lz):
		self.n_x   = nx
		self.n_y   = ny
		self.n_z   = nz
		self.L_x   = Lx
		self.L_y   = Ly
		self.L_z   = Lz
		self.d_t   = dt
 
		self.E_n_x   = np.zeros((nx,ny,nz),dtype=complex)
		self.E_n_y   = np.zeros((nx,ny,nz),dtype=complex)
		self.E_n_z   = np.zeros((nx,ny,nz),dtype=complex)
	
		self.E_old_x = np.zeros((nx,ny,nz),dtype=complex)
		self.E_old_y = np.zeros((nx,ny,nz),dtype=complex)
		self.E_old_z = np.zeros((nx,ny,nz),dtype=complex)

		self.B_n_x = np.zeros((nx,ny,nz),dtype=complex)
		self.B_n_y = np.zeros((nx,ny,nz),dtype=complex)
		self.B_n_z = np.zeros((nx,ny,nz),dtype=complex)

		self.B_old_x = np.zeros((nx,ny,nz),dtype=complex)
		self.B_old_y = np.zeros((nx,ny,nz),dtype=complex)
		self.B_old_z = np.zeros((nx,ny,nz),dtype=complex)

	        self.mesh_x  = np.zeros((nx,ny,nz))
		self.mesh_y  = np.zeros((nx,ny,nz))
		self.mesh_z  = np.zeros((nx,ny,nz))

		self.d_x   = (Lx-0)/nx
	   	self.d_y   = (Ly-0)/ny
		self.d_z   = (Lz-0)/nz

		for i in range (nx):
			for j in range (ny):
				for k in range (nz):
					self.mesh_x[i,j,k]=i*self.d_x
					self.mesh_y[i,j,k]=j*self.d_y
					self.mesh_z[i,j,k]=k*self.d_z		


	def Time_Integration_PSM(self):
		print "test %d " %self.n_x		

