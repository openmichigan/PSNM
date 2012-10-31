#!/usr/bin/env python
"""
Solving Heat Equation using pseudospectral methods with Backwards Euler:
u_t= \alpha*u_xx
BC = u(0)=0 and u(2*pi)=0 (Periodic)
IC=sin(x)
"""

import math
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Grid
N = 64; h = 2*math.pi/N; x = [h*i for i in xrange(1,N+1)]

# Initial conditions
v = [math.sin(y) for y in x]
alpha = 0.5
t = 0     
dt = .001 #Timestep size

# (ik)^2 Vector
I = complex(0,1)
k = numpy.array([I*n for n in range(0,N/2) + [0] + range(-N/2+1,0)])
k2=k**2;

# Setting up Plot
tmax = 5.0; tplot = 0.1
plotgap= int(round(tplot/dt))
nplots = int(round(tmax/tplot))
data = numpy.zeros((nplots+1,N))
data[0,:] = v
tdata = [t]

for i in xrange(nplots):
    v_hat = numpy.fft.fft(v)  # convert to fourier space
    for n in xrange(plotgap):
        v_hat = v_hat / (1-dt*alpha*k2) # backward Euler timestepping

    v = numpy.fft.ifft(v_hat)   # convert back to real space
    data[i+1,:] = numpy.real(v)   # records data

    t = t+plotgap*dt    # records real time
    tdata.append(t)

# Plot using mesh
xx,tt = (numpy.mat(A) for A in (numpy.meshgrid(x,tdata)))

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, tt, data,rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.xlabel('x')
plt.ylabel('t')
plt.show()






