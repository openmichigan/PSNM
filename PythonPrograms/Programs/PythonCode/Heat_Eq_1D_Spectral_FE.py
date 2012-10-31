#!/usr/bin/env python
"""
Solving Heat Equation using pseudo-spectral and Forward Euler
u_t= \alpha*u_xx
BC= u(0)=0, u(2*pi)=0
IC=sin(x)
"""

import math
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Grid
N = 64                     # Number of steps
h = 2*math.pi/N                 # step size
x = h*numpy.arange(0,N)    # discretize x-direction

alpha = 0.5                # Thermal Diffusivity constant
t = 0
dt = .001

# Initial conditions 
v = numpy.sin(x)
I = complex(0,1)
k = numpy.array([I*y for y in range(0,N/2) + [0] + range(-N/2+1,0)])
k2=k**2;

# Setting up Plot
tmax = 5; tplot = .1;
plotgap = int(round(tplot/dt))
nplots  = int(round(tmax/tplot))

data = numpy.zeros((nplots+1,N))
data[0,:] = v
tdata = [t]

for i in xrange(nplots):
    v_hat = numpy.fft.fft(v)

    for n in xrange(plotgap):
        v_hat = v_hat+dt*alpha*k2*v_hat # FE timestepping

    v = numpy.real(numpy.fft.ifft(v_hat))   # back to real space
    data[i+1,:] = v

    # real time vector
    t = t+plotgap*dt
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

