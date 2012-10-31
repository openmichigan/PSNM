"""
A program to solve u_t'=u(u-1) using a Strang
splitting method 
"""

import time
import numpy
import matplotlib.pyplot as plt

Nt = 1000		# Number of timesteps
tmax = 1.0		# Maximum time
dt=tmax/Nt      # increment between times
u0 = 0.8        # Initial value
t0 = 0			# Starting time
u = [u0]		# Variables holding the values of iterations
t = [t0]		# Times of discrete solutions

T0 = time.clock()
for i in xrange(Nt):
	c = -1.0/u[i]
	utemp=-1.0/(c+0.5*dt)
	utemp2=utemp*numpy.exp(-dt)
	c = -1.0/utemp2
	unew=-1.0/(c+0.5*dt)
	u.append(unew)
	t.append(t[i]+dt)

T = time.clock() - T0	
uexact = [4.0/(4.0+numpy.exp(tt)) for tt in t]

print
print "Elapsed time is %f" % (T)

plt.figure()
plt.plot(t,u,'r+',t,uexact,'b-.')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.legend(('Numerical Solution', 'Exact solution'))
plt.title('Numerical solution of du/dt=u(u-1)')
plt.show()