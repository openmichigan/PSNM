#!/usr/bin/env python
"""
Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
u_t= epsilon(u_{xx}+u_{yy}) + u - u^3
where u-u^3 is treated explicitly and u_{xx} and u_{yy} is treated implicitly
BC = Periodic
IC=v=sin(2*pi*x)+0.5*cos(4*pi*y)
"""

import math
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import time
 
plt.ion()

# Setup the grid
N = 64; h = 1.0/N; 
x = [h*i for i in xrange(1,N+1)]
y = [h*i for i in xrange(1,N+1)]
dt = 0.05
xx,yy = (numpy.mat(A) for A in (numpy.meshgrid(x,y)))

# Initial Conditions
u = numpy.array(numpy.sin(2*math.pi*xx) + 0.5*numpy.cos(4*math.pi*yy), dtype=float)

epsilon = 0.01

# (ik) and (ik)^2 vectors in x and y direction
I = complex(0,1)
k_x = numpy.array([I*n for n in range(0,N/2) + [0] + range(-N/2+1,0)])
k_y = k_x

kxx = numpy.zeros((N,N), dtype=complex)
kyy = numpy.zeros((N,N), dtype=complex)
for i in xrange(N):
    for j in xrange(N):
        kxx[i,j] = k_x[i]**2
        kyy[i,j] = k_y[j]**2
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(xx, yy, u,rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

v_hat = numpy.zeros((N,N), dtype=complex)
v_hat = numpy.fft.fft2(u)

for n in xrange(100):
    # calculate nonlinear term in real space
    v_nl = numpy.array(u**3, dtype=complex)        
	# FFT for nonlinear and linear terms
    v_nl = numpy.fft.fft2(v_nl)
    v_hat = (v_hat*(1+1/dt) - v_nl)
    v_hat=v_hat/(1/dt - (kxx+kyy)*epsilon)  # Implicit/Explicit timestepping
    u = numpy.real(numpy.fft.ifft2(v_hat))
    # Remove old plot before drawing
    ax.collections.remove(surf)
    surf = ax.plot_surface(xx, yy, u,rstride=1, cstride=1, cmap=cm.jet, 
         linewidth=0, antialiased=False)
    plt.draw()
plt.show()    
