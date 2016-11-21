import numpy as np
import scipy as sp
import math as math
import cmath as cmath
import matplotlib.pyplot as plt
#resolution de dE/dt=c^2d^2E/dx^2
c=3*math.pow(10,4)
c=1.
#definition du domaine
Lx=1.
nx=128
dx=Lx/nx
dt=0.005
nsteps=100
alpha=c*dt/dx
print "cfl = %f "%alpha
X=np.arange(nx)*Lx*dx
mode=2
w=2*math.pi*c*mode/Lx
k=w/c
print "k=  %f"% k 
E_old=np.zeros(nx,dtype=complex)
E_oldint=np.zeros(nx,dtype=complex)
for i in range(nx):
	E_old[i]=cmath.sin(-2*math.pi*X[i])
plt.plot(X,E_old)
plt.show()
E_n=0*E_old
E_ntilde=0*E_old
E_old_tilde=np.fft.fft(E_old)
x=range(nx/2)
rx=range(-nx/2+1,0,1)
x.append(0)
x.extend(rx)
K_mesh=2*math.pi/Lx*np.asarray(x)

for i in range(nsteps):
	print "time = %f i = %d " %(dt*i,i)
	E_ntilde=E_old_tildedt-K_mesh*K_mesh*E_old_tilde
	E_n=np.fft.ifft(E_ntilde)
	E_old_tilde=E_ntilde
	if (i%5 ==0):
		plt.plot(E_n)
		plt.show()




