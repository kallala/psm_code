import numpy as np
import scipy as sp
import math  as math
import cmath as cmath


def prod_vec(k,v): #calcule le produit scalaire de deux vecteurs
	r=np.zeros(3,dtype=complex)
	r[0]=k[1]*v[2]-k[2]*v[1]
	r[1]=k[2]*v[0]-k[0]*v[1]
	r[2]=k[0]*v[1]-k[1]*v[0]
	return(r)

def prod_vec_champs(K,E,i,j,k): #calcule le  produit scalaire sur un point du maillage entre le maillage spectral et le champs E
	k_=np.zeros(3)
	e=k_
	for ii in range(3):
		k_[ii]=K[ii][j,i,k]
		e[ii]=E[i,j,k,ii]
	return(prod_vec(k,v))
def prod_vec_champs_totaux(K,E):   #calcule le produit scalaire sur tout le maillage entre E et le maillage spectral
	Etilde=E
	nx=E.shape[0]
	ny=E.shape[1]
	nz=E.shape[2]
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				Etilde[i,j,k,:]=prod_vec_champs(K,E,i,j,k)
	return(Etilde)	
def fft_champs_3d (E):
	Etilde = E
 	nx     = E.shape[0]
	ny     = E.shape[1]
	nz     = E.shape[2]
	for d in range(3):
		for i in range(nx):
			for j in range(ny):
				Etilde[i,j,:,d] = np.fft.fft(E[i,j,:,d])
		for j in range(ny):
			for k in range(nz):
				Etilde[:,j,k,d] = np.fft.fft(Etilde[:,j,k,d])		

			for i in range(nx):
				Etilde[i,:,k,d] = np.fft.fft(Etilde[i,:,k,d])
	return(Etilde)	
		
	

def ifft_champs_3d(Etilde):
	E  = Etilde
        nx = Etilde.shape[0]
        ny = Etilde.shape[1]
        nz = Etilde.shape[2]
        for d in range(3):
                for i in range(nx):
                        for j in range(ny):
                                E[i,j,:,d] = np.fft.ifft(Etilde[i,j,:,d])
                for j in range(ny):
                        for k in range(nz):
                                E[:,j,k,d] = np.fft.ifft(E[:,j,k,d])
                for k in range(nz):
                        for i in range(nx):
                                E[i,:,k,d] = np.fft.ifft(E[i,:,k,d])
        return(E)




















