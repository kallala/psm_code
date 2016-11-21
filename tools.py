import numpy as np
import scipy as sp
import math  as math
import cmath as cmath
import matplotlib.pyplot as plt


def prod_vec(k,v): #calcule le produit scalaire de deux vecteurs
	r=np.zeros(3,dtype=complex)
	r[0]=k[1]*v[2]-k[2]*v[1]
	r[1]=-(k[0]*v[2]-k[2]*v[0])
	r[2]=k[0]*v[1]-k[1]*v[0]
	return(r)

def prod_vec_champs(K,E,i,j,k): #calcule le  produit scalaire sur un point du maillage entre le maillage spectral et le champs E
	k_=np.zeros(3,dtype=complex)
  	e=np.zeros(3,dtype=complex)	
	for ii in range(3):
		k_[ii]=K[ii,i,j,k]
		e[ii]=E[ii,i,j,k]
	v=prod_vec(k_,e)
	return(v)
def prod_vec_champs_totaux(K,E):   #calcule le produit scalaire sur tout le maillage entre E et le maillage spectral
	Etilde=np.copy(E)
	nx=E.shape[1]
	ny=E.shape[2]
	nz=E.shape[3]
	Etilde[0,:,:,:]=K[1,:,:,:]*E[2,:,:,:]-E[1,:,:,:]*K[2,:,:,:]
	Etilde[1,:,:,:]=-(K[0,:,:,:]*E[2,:,:,]-E[0,:,:,:]*K[2,:,:,:])
	Etilde[2,:,:,:]=K[0,:,:,:]*E[1,:,:,:]-E[0,:,:,:]*K[1,:,:,:]
	#for i in range(nx):
	#	for j in range(ny):
	#		for k in range(nz):
	#			Etilde[:,i,j,k]=prod_vec_champs(K,E,i,j,k)
	return(Etilde)	

def ifft_champs_3d(Etilde):
	E  = np.copy(Etilde)
        nx = Etilde.shape[1]
        ny = Etilde.shape[2]
        nz = Etilde.shape[3]
        # for d in range(3):
        #        for i in range(nx):
        #                for k in range(nz):
        #                        E[d,i,:,k] = np.fft.ifft(E[d,i,:,k])
        #        for j in range(ny):
        #                for k in range(nz):
        #                        E[d,:,j,k] = np.fft.ifft(E[d,:,j,k])
        #        for j in range(ny):
        #                for i in range(nx):
        #                        E[d,i,j,:] = np.fft.ifft(E[d,i,j,:])
        
	for d in range(3):
		E[d,:,:,:]=np.fft.ifftn(Etilde[d,:,:,:])
	
	return(E)
def egalite(E):
	A=E
	return(A)
def fft_champs_3d(E):
        Etilde=np.copy(E)
	nx = Etilde.shape[0]
        ny = Etilde.shape[1]
        nz = Etilde.shape[2]
        #for d in range(3):
        #        for i in range(nx):
        #                for j in range(ny):
        #                        Etilde[d,i,j,:] = np.fft.fft(Etilde[d,i,j,:])
        #        for j in range(ny):
        #                for k in range(nz):
        #                        Etilde[d,:,j,k] = np.fft.fft(Etilde[d,:,j,k])
        #        for k in range(nz):
        #                for i in range(nx):
        #                        Etilde[d,i,:,k] = np.fft.fft(Etilde[d,i,:,k])
       
	for d in range(3):
		Etilde[d,:,:,:]=np.fft.fftn(E[d,:,:,:])
	return(Etilde)
