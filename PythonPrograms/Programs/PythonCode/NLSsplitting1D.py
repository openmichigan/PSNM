"""
A program to solve the 1D Nonlinear Schrodinger equation using a
second order splitting method. The numerical solution is compared
to an exact soliton solution. The exact equation solved is
iu_t+u_{xx}+2|u|^2u=0

More information on visualization can be found on the Mayavi
website, in particular:
http://github.enthought.com/mayavi/mayavi/mlab.html
which was last checked on 6 April 2012

"""

import math
import numpy
import matplotlib.pyplot as plt
import time

plt.ion()

# Grid
Lx=16.0 	    # Period 2*pi*Lx
Nx=8192 	    # Number of harmonics
Nt=1000 	    # Number of time slices
tmax=1.0    # Maximum time
dt=tmax/Nt  # time step
plotgap=10  # time steps between plots
Es= -1.0     # focusing (+1) or defocusing (-1) parameter
numplots=Nt/plotgap  # number of plots to make

x = [i*2.0*math.pi*(Lx/Nx) for i in xrange(-Nx/2,1+Nx/2)]
k_x = (1.0/Lx)*numpy.array([complex(0,1)*n for n in range(0,Nx/2) \
+ [0] + range(-Nx/2+1,0)])

k2xm=numpy.zeros((Nx), dtype=float)
xx=numpy.zeros((Nx), dtype=float)

for i in xrange(Nx):
   k2xm[i] = numpy.real(k_x[i]**2)
   xx[i]=x[i]
        

# allocate arrays
usquared=numpy.zeros((Nx), dtype=float)
pot=numpy.zeros((Nx), dtype=float)
u=numpy.zeros((Nx), dtype=complex)
uexact=numpy.zeros((Nx), dtype=complex)
una=numpy.zeros((Nx), dtype=complex)
unb=numpy.zeros((Nx), dtype=complex)
v=numpy.zeros((Nx), dtype=complex)
vna=numpy.zeros((Nx), dtype=complex)
vnb=numpy.zeros((Nx), dtype=complex)
mass=numpy.zeros((Nx), dtype=complex)
test=numpy.zeros((numplots-1),dtype=float)
tdata=numpy.zeros((numplots-1), dtype=float)

t=0.0
u=4.0*numpy.exp(complex(0,1.0)*t)*\
   (numpy.cosh(3.0*xx)+3.0*numpy.exp(8.0*complex(0,1.0)*t)*numpy.cosh(xx))\
   /(numpy.cosh(4*xx)+4.0*numpy.cosh(2.0*xx)+3.0*numpy.cos(8.0*t))
uexact=u
v=numpy.fft.fftn(u)
usquared=abs(u)**2
fig =plt.figure()
ax = fig.add_subplot(311)
ax.plot(xx,numpy.real(u),'b-')
plt.xlabel('x')
plt.ylabel('real u')
ax = fig.add_subplot(312)
ax.plot(xx,numpy.imag(u),'b-')
plt.xlabel('x')
plt.ylabel('imaginary u')
ax = fig.add_subplot(313)
ax.plot(xx,abs(u-uexact),'b-')
plt.xlabel('x')
plt.ylabel('error')
plt.show()


# initial mass
usquared=abs(u)**2
mass=numpy.fft.fftn(usquared)
ma=numpy.real(mass[0])
maO=ma
tdata[0]=t
plotnum=0
#solve pde and plot results
for nt in xrange(numplots-1):
    for n in xrange(plotgap):
        vna=v*numpy.exp(complex(0,0.5)*dt*k2xm)
        una=numpy.fft.ifftn(vna)
        usquared=2.0*abs(una)**2
        pot=Es*usquared
        unb=una*numpy.exp(complex(0,-1)*dt*pot)
        vnb=numpy.fft.fftn(unb)
        v=vnb*numpy.exp(complex(0,0.5)*dt*k2xm)
        u=numpy.fft.ifftn(v)
        t+=dt
    plotnum+=1
    usquared=abs(u)**2
    uexact = 4.0*numpy.exp(complex(0,1.0)*t)*\
      (numpy.cosh(3.0*xx)+3.0*numpy.exp(8.0*complex(0,1.0)*t)*numpy.cosh(xx))\
      /(numpy.cosh(4*xx)+4.0*numpy.cosh(2.0*xx)+3.0*numpy.cos(8.0*t))
    ax = fig.add_subplot(311)
    plt.cla()
    ax.plot(xx,numpy.real(u),'b-')
    plt.title(t)
    plt.xlabel('x')
    plt.ylabel('real u')
    ax = fig.add_subplot(312)
    plt.cla()
    ax.plot(xx,numpy.imag(u),'b-')
    plt.xlabel('x')
    plt.ylabel('imaginary u')
    ax = fig.add_subplot(313)
    plt.cla()
    ax.plot(xx,abs(u-uexact),'b-')
    plt.xlabel('x')
    plt.ylabel('error')    
    plt.draw()
    mass=numpy.fft.fftn(usquared)
    ma=numpy.real(mass[0])
    test[plotnum-1]=numpy.log(abs(1-ma/maO))
    print(test[plotnum-1])
    tdata[plotnum-1]=t
       
plt.ioff()
plt.show()