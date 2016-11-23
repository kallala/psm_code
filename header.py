from tools import *

class Simulation :  # class simulation solveur EM 1D pseudo-specral



	def   __init__(self,dt,nx,ny,nz,Lx,Ly,Lz,c):
		self.c = c
		self.n_x     = nx
		self.n_y     = ny
		self.n_z     = nz
		self.L_x     = Lx
		self.L_y     = Ly
		self.L_z     = Lz
		self.d_t     = dt
 
		self.E_n     = np.zeros((3,nx,ny,nz),dtype=complex)
		self.Etilde_n= np.zeros((3,nx,ny,nz),dtype=complex)

		self.E_old   = np.zeros((3,nx,ny,nz),dtype=complex)
		self.Etilde_old = np.zeros((3,nx,ny,nz),dtype=complex)

		self.B_n     = np.zeros((3,nx,ny,nz),dtype=complex)
		self.Btilde_n= np.zeros((3,nx,ny,nz),dtype=complex)
	
		self.B_old   = np.zeros((3,nx,ny,nz),dtype=complex)
		self.Btilde_old = np.zeros((3,nx,ny,nz),dtype=complex)	
		self.d_x     = (Lx-0)/nx
	   	self.d_y     = (Ly-0)/ny
		self.d_z     = (Lz-0)/nz
		#verification de la cfl	
		cfl=self.c*dt/(math.sqrt(self.d_x*self.d_x+self.d_y*self.d_y+self.d_z*self.d_z))
		print " cfl  =   %f  " %(cfl)	
	
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

		[self.mesh_x,self.mesh_y,self.mesh_z]=np.mgrid[0:self.L_x:self.d_x,0:self.L_y:self.d_y,0:self.L_y:self.d_z]
	
		
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
	        self.mesh_K=np.meshgrid(X,Y,Z)
		self.mesh_KK=np.zeros((3,nx,ny,nz))
		for i in range(nx):
			for j in range(ny):
				for k in range(nz):
					for d in range(3):
						self.mesh_KK[d,i,j,k]=self.mesh_K[d][j,i,k]
	#	def fft_champs_3d (self,E):
	#		nx=self.n_x
	#		ny=self.n_y
	#		nz=self.n_z
	#		Etilde=E
	#		for d in range(3):
	#			for i in range(nx):
	#                       	for j in range(ny):
	#                      		Z=Etilde[d,i,j,:]
	#					Etilde[i,j,:,d] = np.fft.fft(Z)
        		      #	for j in range(ny):
                #       	 	for k in range(nz):
                #                	X=Etilde[:,j,k,d]
                #                	Etilde[:,j,k,d] = np.fft.fft(X)         
                #	for k in range(nz):
                #        	for i in range(nx):
                #                	Y=Etilde[i,:,k,d]
		#              		Etilde[i,:,k,d] = np.fft.fft(Y)
	#		return (Etilde)	
	def Initial_Conditions(self,E_0_x,E_0_y,E_0_z,mod,phi_x,phi_y,phi_z):
		c=self.c
		w=2*math.pi*c*mod/self.L_x
		vecteur_onde=w/c
	#	for i in range(self.n_x):
	#		x=self.mesh[0][0,i,0]
	#		print " x = %f " %x
	#		E1=E_0_y*cmath.exp(complex(0,-vecteur_onde*x+phi_y))
	#		print " E 1 = %f + i %f " %( E1.real ,E1.imag)
	#		E2=E_0_z*cmath.exp(complex(0,-vecteur_onde*x+phi_z))
	#		print " E 2 = %f + i %f "%( E2.real,E2.imag)
	#		for j in range(self.n_y):
	#			for k in range( self.n_z):
	#				self.E_old[i,j,k,1]=E1
	#				self.E_old[i,j,k,0]=0
	#				self.E_old[i,j,k,2]=E2
		nx=self.n_x
		ny=self.n_y
		nz=self.n_z
		B_0_y=-E_0_z/c	
		B_0_z=-E_0_y/c
		for i in range(nx):
			for j in range(ny):
				for k in range(nz):
					self.E_old[1,i,j,k]=E_0_y*cmath.exp(complex(0,-vecteur_onde*(self.mesh_y[i,j,k]+self.mesh_x[i,j,k])))
					self.E_old[2,i,j,k]=E_0_z*cmath.exp(complex(0,-vecteur_onde*(self.mesh_x[i,j,k]+self.mesh_y[i,j,k])))
					self.B_old[0,i,j,k]=0
					self.B_old[1,i,j,k]=B_0_y*cmath.exp(complex(0,-vecteur_onde*(self.mesh_y[i,j,k]+self.mesh_x[i,j,k])))
                  		        self.B_old[2,i,j,k]=B_0_z*cmath.exp(complex(0,-vecteur_onde*(self.mesh_y[i,j,k]+self.mesh_x[i,j,k])))
					self.B_old[0,i,j,k]=np.copy(-self.B_old[1,i,j,k])
					self.E_old[0,i,j,k]=np.copy(-self.E_old[1,i,j,k])
		self.Etilde_old=fft_champs_3d(self.E_old)
		self.Btilde_old=fft_champs_3d(self.B_old)
	



	def Time_Integration_PSM(self,nsteps):
		j=complex(0,1)
		c=self.c
		dt=self.d_t
		for i in range(nsteps):	
			if i%100==0:
				print "i= %d " %i
			rotB=j*prod_vec_champs_totaux(self.mesh_KK,self.Btilde_old) #rot B sur le domaine spectral
			self.Etilde_n=self.Etilde_old+dt*c*c*rotB
			rotE=j*prod_vec_champs_totaux(self.mesh_KK,self.Etilde_old)
			self.Btilde_n=-dt*rotE+self.Btilde_old
			self.Btilde_old=np.copy(self.Btilde_n)
			self.Etilde_old=np.copy(self.Etilde_n)
			















