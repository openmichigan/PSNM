#!/usr/bin/env python
"""
A program to solve y'=y^2 using the backward Euler
method and fixed point iteration
This is not optimized and is very simple
"""

import time
import matplotlib.pyplot as plt

N = 1000		# Number of timesteps
tmax = 0.99		# Maximum time
y0 = 1
t0 = 0			# Initial value
tol = pow(0.1,10)	# Tolerance for fixed point iterations
h = tmax/N		# Time step

y = [y0]		# Variables holding the values of iterations
t = [t0]		# Times of discrete solutions



T0 = time.clock()
for i in xrange(N):
	yold = y[i]
	ynew = y[i]
	err = 1
	while err > tol:
		ynew = h*pow(yold,2)+y[i]
		err = abs(ynew-yold)
		yold = ynew
	y.append(ynew)
	t.append(t[i]+h)

T = time.clock() - T0	
yexact = [1.0/(1.0-x) for x in t]

print
print "Exact value:                 y(%d)=%f" % (tmax, 1/(1-tmax))
print "Value given by aproximation: y(%d)=%f" % (tmax, y[-1])
maxerr=(max([abs(y[i] - yexact[i]) for i in xrange(len(y))]))
print "Maximum error:               %f" % maxerr   
print "Elapsed time is %f" % (T)

plt.figure()
plt.plot(t,y,'r+',t,yexact,'b-.')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.legend(('Backward Euler', 'Exact solution'))
plt.title('Numerical solution of dy/dt=y^2')
plt.show()



