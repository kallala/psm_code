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
	def fft_champs_3d (self,E):
		nx=self.n_x
		ny=self.n_y
		nz=self.n_z
		Etilde=E
		for d in range(3):
			for i in range(nx):
                        	for j in range(ny):
                        		Z=E[i,j,:,d]
					Etilde[i,j,:,d] = np.fft.fft(Z)
                	for j in range(ny):
                        	for k in range(nz):
                                	X=Etilde[:,j,k,d]
                                	Etilde[:,j,k,d] = np.fft.fft(X)         
                	for k in range(nz):
                        	for i in range(nx):
                                	Y=Etilde[i,:,k,d]
		              		Etilde[i,:,k,d] = np.fft.fft(Y)
		return (Etilde)	
	def Initial_Conditions(self,E_0_x,E_0_y,E_0_z,mod,phi_x,phi_y,phi_z):
		c=3.*math.pow(10,8)
		w=2*math.pi*c*mod/self.L_x
		vecteur_onde=w/c
		for i in range(self.n_x):
			x=self.mesh[0][0,i,0]
			print " x = %f " %x
			E1=E_0_y*cmath.exp(complex(0,-vecteur_onde*x+phi_y))
			print " E 1 = %f + i %f " %( E1.real ,E1.imag)
			E2=E_0_z*cmath.exp(complex(0,-vecteur_onde*x+phi_z))
			print " E 2 = %f + i %f "%( E2.real,E2.imag)
			for j in range(self.n_y):
				for k in range( self.n_z):
					self.E_old[i,j,k,1]=E1
					self.E_old[i,j,k,0]=0
					self.E_old[i,j,k,2]=E2
		
		Etilde_old=elf.ft_champs_3d(self.E_old)
		Btilde_old=self.fft_champs_3d(self.B_old)
	



	def Time_Integration_PSM(self):
		print "test %d " %self.n_x		
