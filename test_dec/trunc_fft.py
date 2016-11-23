import numpy as np
import scipy as sp
import math  as math
import cmath as cmath
import matplotlib.pyplot as plt

nx=(8)
a=np.zeros(nx)
for i in range(nx):
	x=i*1./nx
	a[i]=math.exp(-x*x/.5)

plt.plot(a)
plt.grid(True)
plt.show()

