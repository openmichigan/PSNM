"""
A program to solve the 1D Klein Gordon equation using a
second order semi-explicit method. The numerical solution is 
compared to an exact solution

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
Lx=64.0 	 # Period 2*pi*Lx
Nx=4096 	 # Number of harmonics
Nt=500 	     # Number of time slices
tmax=5.0     # Maximum time
c=0.5		 # Wave speed
dt=tmax/Nt   # time step
plotgap=10   # time steps between plots
Es= 1.0      # focusing (+1) or defocusing (-1) parameter
numplots=Nt/plotgap  # number of plots to make
tol=0.1**12  # tolerance for fixed point iterations

x = [i*2.0*math.pi*(Lx/Nx) for i in xrange(-Nx/2,1+Nx/2)]
k_x = (1.0/Lx)*numpy.array([complex(0,1)*n for n in range(0,Nx/2) \
+ [0] + range(-Nx/2+1,0)])

kxm=numpy.zeros((Nx), dtype=complex)
xx=numpy.zeros((Nx), dtype=float)

for i in xrange(Nx):
        kxm[i] = k_x[i]
        xx[i] = x[i]
        
# allocate arrays
unew=numpy.zeros((Nx), dtype=float)
u=numpy.zeros((Nx), dtype=float)
utemp=numpy.zeros((Nx), dtype=float)
uexact=numpy.zeros((Nx), dtype=float)
uold=numpy.zeros((Nx), dtype=float)
vnew=numpy.zeros((Nx), dtype=complex)
v=numpy.zeros((Nx), dtype=complex)
vold=numpy.zeros((Nx), dtype=complex)
ux=numpy.zeros((Nx), dtype=float)
vx=numpy.zeros((Nx), dtype=complex)
Kineticenergy=numpy.zeros((Nx), dtype=complex)
Potentialenergy=numpy.zeros((Nx), dtype=complex)
Strainenergy=numpy.zeros((Nx), dtype=complex)
EnKin=numpy.zeros((numplots), dtype=float)
EnPot=numpy.zeros((numplots), dtype=float)
EnStr=numpy.zeros((numplots), dtype=float)
En=numpy.zeros((numplots), dtype=float)
Enchange=numpy.zeros((numplots-1),dtype=float)
tdata=numpy.zeros((numplots), dtype=float)
nonlin=numpy.zeros((Nx), dtype=float)
nonlinhat=numpy.zeros((Nx), dtype=complex)

t=0.0
u=numpy.sqrt(2)/(numpy.cosh((xx-c*t)/numpy.sqrt(1.0-c**2)))
uexact=numpy.sqrt(2)/(numpy.cosh((xx-c*t)/numpy.sqrt(1.0-c**2)))
uold=numpy.sqrt(2)/(numpy.cosh((xx+c*dt)/numpy.sqrt(1.0-c**2)))
v=numpy.fft.fftn(u)
vold=numpy.fft.fftn(uold)
fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot(xx,u,'b-')
plt.xlabel('x')
plt.ylabel('u')
ax=fig.add_subplot(212)
ax.plot(xx,abs(u-uexact),'b-')
plt.xlabel('x')
plt.ylabel('error')
plt.show()
# initial energy
vx=0.5*kxm*(v+vold)
ux=numpy.real(numpy.fft.ifftn(vx))
Kineticenergy=0.5*((u-uold)/dt)**2
Strainenergy=0.5*(ux)**2
Potentialenergy=0.5*(0.5*(u+uold))**2 - Es*0.25*(0.5*(u+uold))**4
Kineticenergy=numpy.fft.fftn(Kineticenergy)
Strainenergy=numpy.fft.fftn(Strainenergy)
Potentialenergy=numpy.fft.fftn(Potentialenergy)
EnKin[0]=numpy.real(Kineticenergy[0])
EnPot[0]=numpy.real(Potentialenergy[0])
EnStr[0]=numpy.real(Strainenergy[0])
En[0]=EnStr[0]+EnPot[0]+EnKin[0]
EnO=En[0]
tdata[0]=t
plotnum=0
#solve pde and plot results
for nt in xrange(numplots-1):
    for n in xrange(plotgap):
        nonlin=(u**2+uold**2)*(u+uold)/4.0
        nonlinhat=numpy.fft.fftn(nonlin)
        chg=1
        unew=u
        while (chg>tol):
            utemp=unew
            vnew=( (0.25*(kxm**2  - 1)*(2*v+vold)\
                 +(2*v-vold)/(dt*dt) +Es*nonlinhat)\
                 /(1/(dt*dt) - (kxm**2  -1)*0.25 ) )
            unew=numpy.real(numpy.fft.ifftn(vnew))
            nonlin=(unew**2+uold**2)*(unew+uold)/4.0
            nonlinhat=numpy.fft.fftn(nonlin)
            chg=numpy.max(abs(unew-utemp))
        t+=dt
        # update old terms
        vold=v
        v=vnew
        uold=u
        u=unew
    plotnum+=1
    uexact=numpy.sqrt(2)/(numpy.cosh((xx-c*t)/numpy.sqrt(1.0-c**2)))
    ax = fig.add_subplot(211)
    plt.cla()
    ax.plot(xx,u,'b-')
    plt.title(t)
    plt.xlabel('x')
    plt.ylabel('u')
    ax = fig.add_subplot(212)
    plt.cla()
    ax.plot(xx,abs(u-uexact),'b-')
    plt.xlabel('x')
    plt.ylabel('error')    
    plt.draw()
    vx=0.5*kxm*(v+vold)
    ux=numpy.real(numpy.fft.ifftn(vx))
    Kineticenergy=0.5*((u-uold)/dt)**2
    Strainenergy=0.5*(ux)**2 
    Potentialenergy=0.5*(0.5*(u+uold))**2 - Es*0.25*(0.5*(u+uold))**4
    Kineticenergy=numpy.fft.fftn(Kineticenergy)
    Strainenergy=numpy.fft.fftn(Strainenergy)
    Potentialenergy=numpy.fft.fftn(Potentialenergy)
    EnKin[plotnum]=numpy.real(Kineticenergy[0])
    EnPot[plotnum]=numpy.real(Potentialenergy[0])
    EnStr[plotnum]=numpy.real(Strainenergy[0])
    En[plotnum]=EnStr[plotnum]+EnPot[plotnum]+EnKin[plotnum]
    Enchange[plotnum-1]=numpy.log(abs(1-En[plotnum]/EnO))
    tdata[plotnum]=t

plt.ioff()

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
