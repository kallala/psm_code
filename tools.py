import numpy as np
import scipy as sp
import math  as math
import cmath as cmath

def multiplication_list(n,List):
	for i in range(len(List)):
		List[i]=n*List[i]
	return List

def prod_vec(k,v):
	r=np.zeros(3,dtype=complex)
	r[0]=k[1]*v[2]-k[2]*v[1]
	r[1]=k[2]*v[0]-k[0]*v[1]
	r[2]=k[0]*v[1]-k[1]*v[0]
	return(r)


