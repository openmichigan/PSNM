"""
A program to solve the 2D Klein Gordon equation using a
second order semi-explicit method 

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
Lx=3.0 	    # Period 2*pi*Lx
Ly=3.0 	    # Period 2*pi*Ly
Nx=512 	    # Number of harmonics
Ny=512	    # Number of harmonics
Nt=200 	    # Number of time slices
tmax=5.0   # Maximum time
dt=tmax/Nt  # time step
plotgap=10  # time steps between plots
Es= 1.0    # focusing (+1) or defocusing (-1) parameter
numplots=Nt/plotgap  # number of plots to make

x = [i*2.0*math.pi*(Lx/Nx) for i in xrange(-Nx/2,1+Nx/2)]
y = [i*2.0*math.pi*(Ly/Ny) for i in xrange(-Ny/2,1+Ny/2)]
k_x = (1.0/Lx)*numpy.array([complex(0,1)*n for n in range(0,Nx/2) \
+ [0] + range(-Nx/2+1,0)])
k_y = (1.0/Ly)*numpy.array([complex(0,1)*n for n in range(0,Ny/2) \
+ [0] + range(-Ny/2+1,0)])

kxm=numpy.zeros((Nx,Ny), dtype=complex)
kym=numpy.zeros((Nx,Ny), dtype=complex)
xx=numpy.zeros((Nx,Ny), dtype=float)
yy=numpy.zeros((Nx,Ny), dtype=float)


for i in xrange(Nx):
    for j in xrange(Ny):
        kxm[i,j] = k_x[i]
        kym[i,j] = k_y[j]
        xx[i,j] = x[i]
        yy[i,j] = y[j]
        

# allocate arrays
unew=numpy.zeros((Nx,Ny), dtype=float)
u=numpy.zeros((Nx,Ny), dtype=float)
uold=numpy.zeros((Nx,Ny), dtype=float)
vnew=numpy.zeros((Nx,Ny), dtype=complex)
v=numpy.zeros((Nx,Ny), dtype=complex)
vold=numpy.zeros((Nx,Ny), dtype=complex)
ux=numpy.zeros((Nx,Ny), dtype=float)
uy=numpy.zeros((Nx,Ny), dtype=float)
vx=numpy.zeros((Nx,Ny), dtype=complex)
vy=numpy.zeros((Nx,Ny), dtype=complex)
Kineticenergy=numpy.zeros((Nx,Ny), dtype=complex)
Potentialenergy=numpy.zeros((Nx,Ny), dtype=complex)
Strainenergy=numpy.zeros((Nx,Ny), dtype=complex)
EnKin=numpy.zeros((numplots), dtype=float)
EnPot=numpy.zeros((numplots), dtype=float)
EnStr=numpy.zeros((numplots), dtype=float)
En=numpy.zeros((numplots), dtype=float)
Enchange=numpy.zeros((numplots-1),dtype=float)
tdata=numpy.zeros((numplots), dtype=float)
nonlin=numpy.zeros((Nx,Ny), dtype=float)
nonlinhat=numpy.zeros((Nx,Ny), dtype=complex)

u=0.1*numpy.exp(-(xx**2 + yy**2))*numpy.sin(10*xx+12*yy) 
uold=u
v=numpy.fft.fft2(u)
vold=numpy.fft.fft2(uold)
src = mlab.surf(xx,yy,u,colormap='YlGnBu',warp_scale='auto')
mlab.scalarbar(object=src)
mlab.xlabel('x',object=src)
mlab.ylabel('y',object=src)
mlab.zlabel('u',object=src)
# initial energy
vx=0.5*kxm*(v+vold)
vy=0.5*kym*(v+vold)
ux=numpy.fft.ifft2(vx)
uy=numpy.fft.ifft2(vy)
Kineticenergy=0.5*((u-uold)/dt)**2
Strainenergy=0.5*(ux)**2 + 0.5*(uy)**2 
Potentialenergy=0.5*(0.5*(u+uold))**2 - Es*0.25*(0.5*(u+uold))**4
Kineticenergy=numpy.fft.fft2(Kineticenergy)
Strainenergy=numpy.fft.fft2(Strainenergy)
Potentialenergy=numpy.fft.fft2(Potentialenergy)
EnKin[0]=numpy.real(Kineticenergy[0,0])
EnPot[0]=numpy.real(Potentialenergy[0,0])
EnStr[0]=numpy.real(Strainenergy[0,0])
En[0]=EnStr[0]+EnPot[0]+EnKin[0]
EnO=En[0]
t=0.0
tdata[0]=t
plotnum=0
#solve pde and plot results
for nt in xrange(numplots-1):
    for n in xrange(plotgap):
        nonlin=u**3
        nonlinhat=numpy.fft.fft2(nonlin)
        vnew=( (0.25*(kxm**2 + kym**2 - 1)*(2*v+vold) 
          +(2*v-vold)/(dt*dt) +Es*nonlinhat)/ 
          (1/(dt*dt) - (kxm**2 + kym**2 -1)*0.25 ) )
        unew=numpy.real(numpy.fft.ifft2(vnew))
        t+=dt
        # update old terms
        vold=v
        v=vnew
        uold=u
        u=unew
    plotnum+=1
    src.mlab_source.scalars = unew
    vx=0.5*kxm*(v+vold)
    vy=0.5*kym*(v+vold)
    ux=numpy.fft.ifft2(vx)
    uy=numpy.fft.ifft2(vy)
    Kineticenergy=0.5*((u-uold)/dt)**2
    Strainenergy=0.5*(ux)**2 + 0.5*(uy)**2 
    Potentialenergy=0.5*(0.5*(u+uold))**2 - Es*0.25*(0.5*(u+uold))**4
    Kineticenergy=numpy.fft.fft2(Kineticenergy)
    Strainenergy=numpy.fft.fft2(Strainenergy)
    Potentialenergy=numpy.fft.fft2(Potentialenergy)
    EnKin[plotnum]=numpy.real(Kineticenergy[0,0])
    EnPot[plotnum]=numpy.real(Potentialenergy[0,0])
    EnStr[plotnum]=numpy.real(Strainenergy[0,0])
    En[plotnum]=EnStr[plotnum]+EnPot[plotnum]+EnKin[plotnum]
    Enchange[plotnum-1]=numpy.log(abs(1-En[plotnum]/EnO))
    tdata[plotnum]=t

       
plt.figure()
plt.plot(tdata,En,'r+',tdata,EnKin,'b:',tdata,EnPot,'g-.',tdata,EnStr,'y--')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend(('Total', 'Kinetic','Potential','Strain'))
plt.title('Time Dependence of Energy Components')
plt.show()

plt.figure()
plt.plot(Enchange,'r-')
plt.title('Time Dependence of Change in Total Energy')
plt.show()
