#!/usr/bin/env python
"""
Numerical solution of the 2D incompressible Magnetohydrodynamics equations on a
Square Domain [0,1]x[0,1] using a Fourier pseudo-spectral method
and Implicit Midpoint rule timestepping. 

A scalar field is also advected by the flow
\theta_t+u\theta_x+v\theta_y=Dtheta(\theta_{xx}+\theta_{yy})

Periodic free-slip boundary conditions and Initial conditions:
    u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
    v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)

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


dt=0.0025; t=0.0; tmax=1.0
#nplots=int(tmax/dt)
Re=100
Rem=100
Dtheta=0.001

u=numpy.zeros((N,N), dtype=float)
v=numpy.zeros((N,N), dtype=float)
u_y=numpy.zeros((N,N), dtype=float)
v_x=numpy.zeros((N,N), dtype=float)
omega=numpy.zeros((N,N), dtype=float)
theta=numpy.zeros((N,N), dtype=float)
thetaold=numpy.zeros((N,N), dtype=float)
thetax=numpy.zeros((N,N), dtype=float)
thetay=numpy.zeros((N,N), dtype=float)
thetacheck=numpy.zeros((N,N), dtype=float)
alpha=numpy.zeros((N,N), dtype=float)
alphaold=numpy.zeros((N,N), dtype=float)
bx=numpy.zeros((N,N), dtype=float)
bxold=numpy.zeros((N,N), dtype=float)
by=numpy.zeros((N,N), dtype=float)
byold=numpy.zeros((N,N), dtype=float)


# Initial conditions
for i in range(len(x)):
    for j in range(len(y)):
        u[i][j]=numpy.sin(2*math.pi*x[i])*numpy.cos(2*math.pi*y[j])
        v[i][j]=-numpy.cos(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        u_y[i][j]=-2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        v_x[i][j]=2*math.pi*numpy.sin(2*math.pi*x[i])*numpy.sin(2*math.pi*y[j])
        omega[i][j]=v_x[i][j]-u_y[i][j]
	theta[i][j]=numpy.exp(-2.0*((4*math.pi*(x[i]-0.3))**2+(4*math.pi*(y[j]-0.3))**2))
	alpha[i][j]=numpy.sin(4*math.pi*x[i])*numpy.cos(6*math.pi*y[j])


src = mlab.imshow(xx,yy,theta,colormap='jet')
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
phihat=numpy.zeros((N,N), dtype=complex)
alphahat=numpy.zeros((N,N), dtype=complex)
alphahatold=numpy.zeros((N,N), dtype=complex)
omegahat=numpy.zeros((N,N), dtype=complex)
omegahatold=numpy.zeros((N,N), dtype=complex)
thetahat=numpy.zeros((N,N), dtype=complex)
thetahatold=numpy.zeros((N,N), dtype=complex)
nlhat=numpy.zeros((N,N), dtype=complex)
dpsix=numpy.zeros((N,N), dtype=float)
dpsiy=numpy.zeros((N,N), dtype=float)
dphix=numpy.zeros((N,N), dtype=float)
dphiy=numpy.zeros((N,N), dtype=float)
omegacheck=numpy.zeros((N,N), dtype=float)
omegaold=numpy.zeros((N,N), dtype=float)
temp=numpy.zeros((N,N), dtype=float)
omegahat=numpy.fft.fft2(omega)
thetahat=numpy.fft.fft2(theta)
alphahat=numpy.fft.fft2(alpha)

while (t<=tmax):
    chg=1.0

    # Save old values
    uold=u
    vold=v
    omegaold=omega
    omegacheck=omega
    omegahatold = omegahat
    thetaold=theta
    thetahatold=thetahat
    thetacheck=theta
    bxold=bx
    byold=by
    alphahatold=alphahat
    alphaold=alpha
    alphacheck=alpha
    while(chg>tol):
	# fluid field
        # nolinear {n+1,k}
        nlhat=0.25*numpy.fft.fft2((u+uold)*numpy.fft.ifft2((omegahat+omegahatold)*kx)+\
        (v+vold)*numpy.fft.ifft2((omegahat+omegahatold)*ky)-\
         (bx+bxold)*numpy.fft.ifft2((alphahat+alphahatold)*kx)-\
         (by+byold)*numpy.fft.ifft2((alphahat+alphahatold)*ky))

        # Implicit midpoint rule timestepping
        omegahat=((1/dt + 0.5*(1/Re)*(kxx+kyy))*omegahatold \
        -nlhat) \
        /(1/dt -0.5*(1/Re)*(kxx+kyy))

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
	# magnetic field
        nlhat=0.25*numpy.fft.fft2((u+uold)*numpy.fft.ifft2((alphahat+alphahatold)*kx)+\
        (v+vold)*numpy.fft.ifft2((alphahat+alphahatold)*ky)-\
         (bx+bxold)*numpy.fft.ifft2((omegahat+omegahatold)*kx)-\
         (by+byold)*numpy.fft.ifft2((omegahat+omegahatold)*ky))

        # Implicit midpoint rule timestepping
        alphhat=((1/dt + 0.5*(1/Rem)*(kxx+kyy))*alphahatold \
        -nlhat) \
        /(1/dt -0.5*(1/Rem)*(kxx+kyy))

        phihat=-alphahat/(kxx+kyy)
        phihat[0][0]=0
        phihat[N/2][N/2]=0
        phihat[N/2][0]=0
        phihat[0][N/2]=0

        dphix = numpy.real(numpy.fft.ifft2(phihat*kx))
        dphiy = numpy.real(numpy.fft.ifft2(phihat*ky))
        bx=dphiy
        by=-1.0*dphix
        
        alpha=numpy.real(numpy.fft.ifft2(alphahat))
	
	# passive scalar
	thetax=0.5*numpy.real(numpy.fft.ifft2(kx*(thetahat+thetahatold)))
	thetay=0.5*numpy.real(numpy.fft.ifft2(ky*(thetahat+thetahatold)))
	thetahat= numpy.fft.fft2(thetaold-dt*0.5*((uold+u)*thetax+(vold+v)*thetay) ) 
	theta=numpy.real(numpy.fft.ifft2( thetahat))

        temp=abs(omega-omegacheck)
        chg=numpy.max(temp)
	temp=abs(theta-thetacheck)
	chg=chg+numpy.max(temp)
        print(chg)
        omegacheck=omega
	thetacheck=theta
    t+=dt
    src.mlab_source.scalars = theta
