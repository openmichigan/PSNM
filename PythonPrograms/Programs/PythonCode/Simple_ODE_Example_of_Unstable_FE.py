#!/usr/bin/env python
"""
A program to demonstrate instability of timestepping methods# 
when the timestep is inappropriately choosen.################
"""

from math import exp
import matplotlib.pyplot as plt
import numpy

#Differential equation: y'(t)=-l*y(t) y(t_0)=y_0
#Initial Condition, y(t_0)=1 where t_0=0

# Definition of the Grid
h = 0.1					# Time step size
t0 = 0					# Initial value
tmax = 4					# Value to be computed y(tmax)
Npoints = int((tmax-t0)/h)	# Number of points in the Grid

t = [t0]

# Initial data
l = 0.1
y0 = 1			# Initial condition y(t0)=y0
y_be = [y0]		# Variables holding the value given by the Backward Euler Iteration
y_fe = [y0]		# Variables holding the value given by the Forward Euler Iteration
y_imr = [y0]		# Variables holding the value given by the Midpoint Rule Iteration

for i in xrange(Npoints):
    y_fe.append(y_be[-1]*(1-l*h))
    y_be.append(y_fe[-1]/(1+l*h))
    y_imr.append(y_imr[-1]*(2-l*h)/(2+l*h))
    t.append(t[-1]+h)


print
print "Exact Value:          y(%d)=%f" % (tmax, exp(-4))
print "Backward Euler Value: y(%d)=%f" % (tmax, y_be[-1])
print "Forward Euler Value:  y(%d)=%f" % (tmax, y_fe[-1])
print "Midpoint Rule Value:  y(%d)=%f" % (tmax, y_imr[-1])

# Exact Solution
tt=numpy.arange(0,tmax,0.001)
exact = numpy.exp(-l*tt)

# Plot
plt.figure()
plt.plot(tt,exact,'r-',t,y_fe,'b:',t,y_be,'g--',t,y_imr,'k-.');
plt.xlabel('time')
plt.ylabel('y')
plt.legend(('Exact','Forward Euler','Backward Euler',
		'Implicit Midpoint Rule'))
plt.show()



