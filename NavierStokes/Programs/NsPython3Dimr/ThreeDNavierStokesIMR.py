#!/usr/bin/env python
"""
A program to solve the 3D Navier Stokes equations using the
implicit midpoint rule 

The program is based on the Orszag-Patterson algorithm as documented on pg. 98 
of C. Canuto, M.Y. Hussaini, A. Quarteroni and T.A. Zhang 
"Spectral Methods: Evolution to Complex Geometries and Applications to Fluid Dynamics" 
Springer (2007)

The exact solution used to check the numerical method is in
A. Shapiro "The use of an exact solution of the Navier-Stokes equations
in a validation test of a three-dimensional nonhydrostatic numerical
model" Monthly Weather Review vol. 121 pp. 2420-2425 (1993)

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
Lx=1.0 	    	# Period 2*pi*Lx
Ly=1.0 	    	# Period 2*pi*Ly
Lz=1.0 	    	# Period 2*pi*Lz
Nx=64 	    	# Number of harmonics
Ny=64	    	# Number of harmonics
Nz=64 	    	# Number of harmonics
Nt=20 	    	# Number of time slices
tmax=0.2   	# Maximum time
dt=tmax/Nt  	# time step
t=0.0	    	# initial time
Rey=1.0      	# Reynolds number
tol=0.1**(10) 	# tolerance for fixed point iterations

x = [i*2.0*math.pi*(Lx/Nx) for i in xrange(-Nx/2,1+Nx/2)]
y = [i*2.0*math.pi*(Ly/Ny) for i in xrange(-Ny/2,1+Ny/2)]
z = [i*2.0*math.pi*(Lz/Nz) for i in xrange(-Nz/2,1+Nz/2)]
k_x = (1.0/Lx)*numpy.array([complex(0,1)*n for n in range(0,Nx/2) \
+ [0] + range(-Nx/2+1,0)])
k_y = (1.0/Ly)*numpy.array([complex(0,1)*n for n in range(0,Ny/2) \
+ [0] + range(-Ny/2+1,0)])
k_z = (1.0/Lz)*numpy.array([complex(0,1)*n for n in range(0,Nz/2) \
+ [0] + range(-Nz/2+1,0)])

kxm=numpy.zeros((Nx,Ny,Nz), dtype=complex)
kym=numpy.zeros((Nx,Ny,Nz), dtype=complex)
kzm=numpy.zeros((Nx,Ny,Nz), dtype=complex)
k2xm=numpy.zeros((Nx,Ny,Nz), dtype=float)
k2ym=numpy.zeros((Nx,Ny,Nz), dtype=float)
k2zm=numpy.zeros((Nx,Ny,Nz), dtype=float)
xx=numpy.zeros((Nx,Ny,Nz), dtype=float)
yy=numpy.zeros((Nx,Ny,Nz), dtype=float)
zz=numpy.zeros((Nx,Ny,Nz), dtype=float)


for i in xrange(Nx):
    for j in xrange(Ny):
        for k in xrange(Nz):
            kxm[i,j,k] = k_x[i]
            kym[i,j,k] = k_y[j]
            kzm[i,j,k] = k_z[k]
            k2xm[i,j,k] = numpy.real(k_x[i]**2)
            k2ym[i,j,k] = numpy.real(k_y[j]**2)
            k2zm[i,j,k] = numpy.real(k_z[k]**2)
            xx[i,j,k] = x[i]
            yy[i,j,k] = y[j]
            zz[i,j,k] = z[k]
        

# allocate arrays
u=numpy.zeros((Nx,Ny,Nz), dtype=float)
uold=numpy.zeros((Nx,Ny,Nz), dtype=float)
v=numpy.zeros((Nx,Ny,Nz), dtype=float)
vold=numpy.zeros((Nx,Ny,Nz), dtype=float)
w=numpy.zeros((Nx,Ny,Nz), dtype=float)
wold=numpy.zeros((Nx,Ny,Nz), dtype=float)

uexact=numpy.zeros((Nx,Ny,Nz), dtype=float)
vexact=numpy.zeros((Nx,Ny,Nz), dtype=float)
wexact=numpy.zeros((Nx,Ny,Nz), dtype=float)

utemp=numpy.zeros((Nx,Ny,Nz), dtype=float)
vtemp=numpy.zeros((Nx,Ny,Nz), dtype=float)
wtemp=numpy.zeros((Nx,Ny,Nz), dtype=float)

omegax=numpy.zeros((Nx,Ny,Nz), dtype=float)
omegay=numpy.zeros((Nx,Ny,Nz), dtype=float)
omegaz=numpy.zeros((Nx,Ny,Nz), dtype=float)
omegatot=numpy.zeros((Nx,Ny,Nz), dtype=float)

ux=numpy.zeros((Nx,Ny,Nz), dtype=float)
uy=numpy.zeros((Nx,Ny,Nz), dtype=float)
uz=numpy.zeros((Nx,Ny,Nz), dtype=float)
vx=numpy.zeros((Nx,Ny,Nz), dtype=float)
vy=numpy.zeros((Nx,Ny,Nz), dtype=float)
vz=numpy.zeros((Nx,Ny,Nz), dtype=float)
wx=numpy.zeros((Nx,Ny,Nz), dtype=float)
wy=numpy.zeros((Nx,Ny,Nz), dtype=float)
wz=numpy.zeros((Nx,Ny,Nz), dtype=float)

uxold=numpy.zeros((Nx,Ny,Nz), dtype=float)
uyold=numpy.zeros((Nx,Ny,Nz), dtype=float)
uzold=numpy.zeros((Nx,Ny,Nz), dtype=float)
vxold=numpy.zeros((Nx,Ny,Nz), dtype=float)
vyold=numpy.zeros((Nx,Ny,Nz), dtype=float)
vzold=numpy.zeros((Nx,Ny,Nz), dtype=float)
wxold=numpy.zeros((Nx,Ny,Nz), dtype=float)
wyold=numpy.zeros((Nx,Ny,Nz), dtype=float)
wzold=numpy.zeros((Nx,Ny,Nz), dtype=float)

nonlinu=numpy.zeros((Nx,Ny,Nz), dtype=float)
nonlinv=numpy.zeros((Nx,Ny,Nz), dtype=float)
nonlinw=numpy.zeros((Nx,Ny,Nz), dtype=float)

uhat=numpy.zeros((Nx,Ny,Nz), dtype=complex)
what=numpy.zeros((Nx,Ny,Nz), dtype=complex)
vhat=numpy.zeros((Nx,Ny,Nz), dtype=complex)
phat=numpy.zeros((Nx,Ny,Nz), dtype=complex)
temphat=numpy.zeros((Nx,Ny,Nz), dtype=complex)

rhsuhatfix=numpy.zeros((Nx,Ny,Nz), dtype=complex)
rhswhatfix=numpy.zeros((Nx,Ny,Nz), dtype=complex)
rhsvhatfix=numpy.zeros((Nx,Ny,Nz), dtype=complex)

nonlinuhat=numpy.zeros((Nx,Ny,Nz), dtype=complex)
nonlinvhat=numpy.zeros((Nx,Ny,Nz), dtype=complex)
nonlinwhat=numpy.zeros((Nx,Ny,Nz), dtype=complex)

tdata=numpy.zeros((Nt), dtype=float)

# initial conditions for Taylor-Green Vortex
#theta=0.0
#u=(2.0/(3.0**0.5))*numpy.sin(theta+2.0*math.pi/3.0)*numpy.sin(xx)*numpy.cos(yy)*numpy.cos(zz)
#v=(2.0/(3.0**0.5))*numpy.sin(theta-2.0*math.pi/3.0)*numpy.cos(xx)*numpy.sin(yy)*numpy.cos(zz)
#w=(2.0/(3.0**0.5))*numpy.sin(theta)*numpy.cos(xx)*numpy.cos(yy)*numpy.sin(zz)

# Exact solution
sl=1
sk=1
sm=1
lamlkm=(sl**2.0+sk**2.0+sm**2.0)**0.5 
u=-0.5*(lamlkm*sl*numpy.cos(sk*xx)*numpy.sin(sl*yy)*numpy.sin(sm*zz) \
+sm*sk*numpy.sin(sk*xx)*numpy.cos(sl*yy)*numpy.cos(sm*zz))*numpy.exp(-t*(lamlkm**2.0)/Rey)
v= 0.5*(lamlkm*sk*numpy.sin(sk*xx)*numpy.cos(sl*yy)*numpy.sin(sm*zz) \
-sm*sl*numpy.cos(sk*xx)*numpy.sin(sl*yy)*numpy.cos(sm*zz))*numpy.exp(-t*(lamlkm**2.0)/Rey)
w= numpy.cos(sk*xx)*numpy.cos(sl*yy)*numpy.sin(sm*zz)*numpy.exp(-t*(lamlkm**2.0)/Rey)

uhat=numpy.fft.fftn(u)
vhat=numpy.fft.fftn(v)
what=numpy.fft.fftn(w)

temphat=kxm*uhat
ux=numpy.real(numpy.fft.ifftn(temphat))
temphat=kym*uhat
uy=numpy.real(numpy.fft.ifftn(temphat))
temphat=kzm*uhat
uz=numpy.real(numpy.fft.ifftn(temphat))

temphat=kxm*vhat
vx=numpy.real(numpy.fft.ifftn(temphat))
temphat=kym*vhat
vy=numpy.real(numpy.fft.ifftn(temphat))
temphat=kzm*vhat
vz=numpy.real(numpy.fft.ifftn(temphat))

temphat=kxm*what
wx=numpy.real(numpy.fft.ifftn(temphat))
temphat=kym*what
wy=numpy.real(numpy.fft.ifftn(temphat))
temphat=kzm*what
wz=numpy.real(numpy.fft.ifftn(temphat))

# Calculate vorticity for plotting

omegax=wy-vz
omegay=uz-wx
omegaz=vx-uy
omegatot=omegax**2.0 + omegay**2.0 + omegaz**2.0


#src=mlab.contour3d(xx,yy,zz,u,colormap='jet',opacity=0.1,contours=4)
src = mlab.pipeline.scalar_field(xx,yy,zz,omegatot,colormap='YlGnBu')
mlab.pipeline.iso_surface(src, contours=[omegatot.min()+0.1*omegatot.ptp(), ], \
   colormap='YlGnBu',opacity=0.85)
mlab.pipeline.iso_surface(src, contours=[omegatot.max()-0.1*omegatot.ptp(), ], \
   colormap='YlGnBu',opacity=1.0)
mlab.pipeline.image_plane_widget(src,plane_orientation='z_axes', \
                            slice_index=Nz/2,colormap='YlGnBu', \
                            opacity=0.01)
mlab.pipeline.image_plane_widget(src,plane_orientation='y_axes', \
                            slice_index=Ny/2,colormap='YlGnBu', \
                            opacity=0.01)
mlab.pipeline.image_plane_widget(src,plane_orientation='x_axes', \
                            slice_index=Nx/2,colormap='YlGnBu', \
                            opacity=0.01)
mlab.scalarbar()
mlab.xlabel('x',object=src)
mlab.ylabel('y',object=src)
mlab.zlabel('z',object=src)



t=0.0
tdata[0]=t
#solve pde and plot results
for n in xrange(Nt):
	uold=u
	uxold=ux
	uyold=uy
	uzold=uz
	vold=v
	vxold=vx
	vyold=vy
	vzold=vz
	wold=w
	wxold=wx
	wyold=wy
	wzold=wz
	rhsuhatfix=(1.0/dt + (0.5/Rey)*(k2xm+k2ym+k2zm))*uhat
	rhsvhatfix=(1.0/dt + (0.5/Rey)*(k2xm+k2ym+k2zm))*vhat
	rhswhatfix=(1.0/dt + (0.5/Rey)*(k2xm+k2ym+k2zm))*what
	chg=1.0
	t=t+dt
	while(chg>tol):
		nonlinu=0.25*((u+uold)*(ux+uxold)+(v+vold)*(uy+uyold)+(w+wold)*(uz+uzold))
		nonlinv=0.25*((u+uold)*(vx+vxold)+(v+vold)*(vy+vyold)+(w+wold)*(vz+vzold))
		nonlinw=0.25*((u+uold)*(wx+wxold)+(v+vold)*(wy+wyold)+(w+wold)*(wz+wzold))
		nonlinuhat=numpy.fft.fftn(nonlinu)
		nonlinvhat=numpy.fft.fftn(nonlinv)
		nonlinwhat=numpy.fft.fftn(nonlinw)
		phat=-1.0*(kxm*nonlinuhat+kym*nonlinvhat+kzm*nonlinwhat)/(k2xm+k2ym+k2zm+0.1**13)
		uhat=(rhsuhatfix-nonlinuhat-kxm*phat)/(1.0/dt - (0.5/Rey)*(k2xm+k2ym+k2zm))
		vhat=(rhsvhatfix-nonlinvhat-kym*phat)/(1.0/dt - (0.5/Rey)*(k2xm+k2ym+k2zm))
		what=(rhswhatfix-nonlinwhat-kzm*phat)/(1.0/dt - (0.5/Rey)*(k2xm+k2ym+k2zm))

		temphat=kxm*uhat
		ux=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kym*uhat
		uy=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kzm*uhat
		uz=numpy.real(numpy.fft.ifftn(temphat))

		temphat=kxm*vhat
		vx=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kym*vhat
		vy=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kzm*vhat
		vz=numpy.real(numpy.fft.ifftn(temphat))

		temphat=kxm*what
		wx=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kym*what
		wy=numpy.real(numpy.fft.ifftn(temphat))
		temphat=kzm*what
		wz=numpy.real(numpy.fft.ifftn(temphat))

		utemp=u
		vtemp=v
		wtemp=w
		u=numpy.real(numpy.fft.ifftn(uhat))
		v=numpy.real(numpy.fft.ifftn(vhat))
		w=numpy.real(numpy.fft.ifftn(what))

		chg=numpy.max(abs(u-utemp))+numpy.max(abs(v-vtemp))+numpy.max(abs(w-wtemp))
		
	# calculate vorticity for plotting
	omegax=wy-vz
	omegay=uz-wx
	omegaz=vx-uy
	omegatot=omegax**2.0 + omegay**2.0 + omegaz**2.0
	src.mlab_source.scalars = omegatot
	tdata[n]=t

       
uexact=-0.5*(lamlkm*sl*numpy.cos(sk*xx)*numpy.sin(sl*yy)*numpy.sin(sm*zz) \
+ sm*sk*numpy.sin(sk*xx)*numpy.cos(sl*yy)*numpy.cos(sm*zz))*numpy.exp(-t*(lamlkm**2.0)/Rey)
vexact= 0.5*(lamlkm*sk*numpy.sin(sk*xx)*numpy.cos(sl*yy)*numpy.sin(sm*zz) \
- sm*sl*numpy.cos(sk*xx)*numpy.sin(sl*yy)*numpy.cos(sm*zz))*numpy.exp(-t*(lamlkm**2.0)/Rey)
wexact= numpy.cos(sk*xx)*numpy.cos(sl*yy)*numpy.sin(sm*zz)*numpy.exp(-t*(lamlkm**2.0)/Rey)

err=numpy.max(abs(u-uexact))+numpy.max(abs(v-vexact))+numpy.max(abs(w-wexact))
print(err)
