#!/usr/bin/env python
"""
Numerical solution of the 2D incompressible Navier-Stokes on a
Square Domain [0,1]x[0,1] unumpy.sing a Fourier pseudo-spectral method
and Crank-Nicolson timestepmath.ping. The numerical solution is compared to
the exact Taylor-Green Vortex solution of the Navier-Stokes equations

Periodic free-slip boundary conditions and Initial conditions:
    u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
    v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
Analytical Solution:
    u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*nu*t)
    v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*nu*t)
"""

import math
import numpy
import matplotlib.pyplot as plt
from mayavi import mlab
import time

# Grid
N=64; h=1.0/N
x = [h*i for i in xrange(1,N+1)]
y = [h*i for i in xrange(1,N+1)]
numpy.savetxt('x.txt',x)

xx=numpy.zeros((N,N), dtype=float)
yy=numpy.zeros((N,N), dtype=float)

for i in xrange(N):
    for j in xrange(N):
        xx[i,j] = x[i]
        yy[i,j] = y[j]


dt=0.0025; t=0.0; tmax=0.10
#nplots=int(tmax/dt)
Rey=1

u=numpy.zeros((N,N), dtype=float)
v=numpy.zeros((N,N), dtype=float)
u_y=numpy.zeros((N,N), dtype=float)
v_x=numpy.zeros((N,N), dtype=float)
omega=numpy.zeros((N,N), dtype=float)
# Initial conditions
for i in range(len(x)):
    for j in range(len(y)):
        u[i][j]=numpy.sin(2*math.pi*x[i])*numpy.cos(2*math.pi*y[j])
        v[i][j]=-numpy.cos(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        u_y[i][j]=-2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        v_x[i][j]=2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        omega[i][j]=v_x[i][j]-u_y[i][j]

src = mlab.imshow(xx,yy,omega,colormap='jet')
mlab.scalarbar(object=src)
mlab.xlabel('x',object=src)
mlab.ylabel('y',object=src)


# Wavenumber
k_x = 2*math.pi*numpy.array([complex(0,1)*n for n in range(0,N/2) \
+ [0] + range(-N/2+1,0)])
k_y=k_x

kx=numpy.zeros((N,N), dtype=complex)
ky=numpy.zeros((N,N), dtype=complex)
kxx=numpy.zeros((N,N), dtype=complex)
kyy=numpy.zeros((N,N), dtype=complex)

for i in xrange(N):
    for j in xrange(N):
        kx[i,j] = k_x[i]
        ky[i,j] = k_y[j]
        kxx[i,j] = k_x[i]**2
        kyy[i,j] = k_y[j]**2

tol=10**(-10)
psihat=numpy.zeros((N,N), dtype=complex)
omegahat=numpy.zeros((N,N), dtype=complex)
omegahatold=numpy.zeros((N,N), dtype=complex)
nlhat=numpy.zeros((N,N), dtype=complex)
nlhatold=numpy.zeros((N,N), dtype=complex)
dpsix=numpy.zeros((N,N), dtype=float)
dpsiy=numpy.zeros((N,N), dtype=float)
omegacheck=numpy.zeros((N,N), dtype=float)
omegaold=numpy.zeros((N,N), dtype=float)
temp=numpy.zeros((N,N), dtype=float)
omegahat=numpy.fft.fft2(omega)
nlhat=numpy.fft.fft2(u*numpy.fft.ifft2(omegahat*kx)+\
v*numpy.fft.ifft2(omegahat*ky))
while (t<=tmax):
    chg=1.0

    # Save old values
    uold=u
    vold=v
    omegaold=omega
    omegacheck=omega
    omegahatold = omegahat
    nlhatold=nlhat
    
    while(chg>tol):
        # nolinear {n+1,k}
        nlhat=numpy.fft.fft2(u*numpy.fft.ifft2(omegahat*kx)+\
        v*numpy.fft.ifft2(omegahat*ky))

        # Crank-Nicolson timestepmath.ping
        omegahat=((1/dt + 0.5*(1/Rey)*(kxx+kyy))*omegahatold \
        -0.5*(nlhatold+nlhat)) \
        /(1/dt -0.5*(1/Rey)*(kxx+kyy))

        psihat=-omegahat/(kxx+kyy)
        psihat[0][0]=0
        psihat[N/2][N/2]=0
        psihat[N/2][0]=0
        psihat[0][N/2]=0

        dpsix = numpy.real(numpy.fft.ifft2(psihat*kx))
        dpsiy = numpy.real(numpy.fft.ifft2(psihat*ky))
        u=dpsiy
        v=-1.0*dpsix
        
        omega=numpy.real(numpy.fft.ifft2(omegahat))
        temp=abs(omega-omegacheck)
        chg=numpy.max(temp)
        print(chg)
        omegacheck=omega
    t+=dt
    src.mlab_source.scalars = omega
 
omegaexact=numpy.zeros((N,N), dtype=float)
for i in range(len(x)):
    for j in range(len(y)):
        uexact_y=-2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*x[j])\
        *numpy.exp(-8*(math.pi**2)*t/Rey)
        vexact_x=2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])\
        *numpy.exp(-8*(math.pi**2)*t/Rey)
        omegaexact[i][j]=vexact_x-uexact_y
numpy.savetxt('Error.txt',abs(omegaexact-omega))




















