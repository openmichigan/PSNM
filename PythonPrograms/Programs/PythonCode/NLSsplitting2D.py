"""
A program to solve the 2D Nonlinear Schrodinger equation using a
second order splitting method 

More information on visualization can be found on the Mayavi
website, in particular:
http://github.enthought.com/mayavi/mayavi/mlab.html
which was last checked on 6 April 2012

"""

import math
import numpy
from mayavi import mlab
import matplotlib.pyplot as plt
import time

# Grid
Lx=4.0 	    # Period 2*pi*Lx
Ly=4.0 	    # Period 2*pi*Ly
Nx=64 	    # Number of harmonics
Ny=64 	    # Number of harmonics
Nt=100 	    # Number of time slices
tmax=1.0    # Maximum time
dt=tmax/Nt  # time step
plotgap=10  # time steps between plots
Es= 1.0     # focusing (+1) or defocusing (-1) parameter
numplots=Nt/plotgap  # number of plots to make

x = [i*2.0*math.pi*(Lx/Nx) for i in xrange(-Nx/2,1+Nx/2)]
y = [i*2.0*math.pi*(Ly/Ny) for i in xrange(-Ny/2,1+Ny/2)]
k_x = (1.0/Lx)*numpy.array([complex(0,1)*n for n in range(0,Nx/2) \
+ [0] + range(-Nx/2+1,0)])
k_y = (1.0/Ly)*numpy.array([complex(0,1)*n for n in range(0,Ny/2) \
+ [0] + range(-Ny/2+1,0)])

k2xm=numpy.zeros((Nx,Ny), dtype=float)
k2ym=numpy.zeros((Nx,Ny), dtype=float)
xx=numpy.zeros((Nx,Ny), dtype=float)
yy=numpy.zeros((Nx,Ny), dtype=float)


for i in xrange(Nx):
    for j in xrange(Ny):
            k2xm[i,j] = numpy.real(k_x[i]**2)
            k2ym[i,j] = numpy.real(k_y[j]**2)
            xx[i,j]=x[i]
            yy[i,j]=y[j]
        

# allocate arrays
usquared=numpy.zeros((Nx,Ny), dtype=float)
pot=numpy.zeros((Nx,Ny), dtype=float)
u=numpy.zeros((Nx,Ny), dtype=complex)
una=numpy.zeros((Nx,Ny), dtype=complex)
unb=numpy.zeros((Nx,Ny), dtype=complex)
v=numpy.zeros((Nx,Ny), dtype=complex)
vna=numpy.zeros((Nx,Ny), dtype=complex)
vnb=numpy.zeros((Nx,Ny), dtype=complex)
mass=numpy.zeros((Nx,Ny), dtype=complex)
test=numpy.zeros((numplots-1),dtype=float)
tdata=numpy.zeros((numplots-1), dtype=float)

u=numpy.exp(-(xx**2 + yy**2 )) 
v=numpy.fft.fftn(u)
usquared=abs(u)**2
src = mlab.surf(xx,yy,usquared,colormap='YlGnBu',warp_scale='auto')
mlab.scalarbar()
mlab.xlabel('x',object=src)
mlab.ylabel('y',object=src)
mlab.zlabel('abs(u)^2',object=src)

# initial mass
usquared=abs(u)**2
mass=numpy.fft.fftn(usquared)
ma=numpy.real(mass[0,0])
print(ma)
maO=ma
t=0.0
tdata[0]=t
plotnum=0
#solve pde and plot results
for nt in xrange(numplots-1):
    for n in xrange(plotgap):
        vna=v*numpy.exp(complex(0,0.5)*dt*(k2xm+k2ym))
        una=numpy.fft.ifftn(vna)
        usquared=abs(una)**2
        pot=Es*usquared
        unb=una*numpy.exp(complex(0,-1)*dt*pot)
        vnb=numpy.fft.fftn(unb)
        v=vnb*numpy.exp(complex(0,0.5)*dt*(k2xm+k2ym) )
        u=numpy.fft.ifftn(v)
        t+=dt
    plotnum+=1
    usquared=abs(u)**2
    src.mlab_source.scalars = usquared
    mass=numpy.fft.fftn(usquared)
    ma=numpy.real(mass[0,0])
    test[plotnum-1]=numpy.log(abs(1-ma/maO))
    print(test[plotnum-1])
    tdata[plotnum-1]=t
       
plt.figure()
plt.plot(tdata,test,'r-')
plt.title('Time Dependence of Change in Mass')
plt.show()
