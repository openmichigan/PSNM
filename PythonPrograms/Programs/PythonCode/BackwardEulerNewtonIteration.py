#!/usr/bin/env python
"""
A program to solve y'=y^2 using the backward Euler
method and Newton iteration
This is not optimized and is very simple
"""

import time
import matplotlib.pyplot as plt

N = 1000            # Number of timesteps
tmax = 0.99         # Maximum value
t0 = 0              # Initial t value
y0 = 1              # Initial value y(t0) = y0
tol = 0.1**10       # Tolerance for fixed point iterations
h = (tmax - t0)/N   # Time step

y = [y0]            # List for discrete solutions
t = [t0]            # List with grid of discrete values of t

T0 = time.clock()       #Start timing

for i in xrange(N):
    yold = y[i]
    ynew = y[i]
    err =1
    while err > tol:
        ynew = yold-(yold-y[i]-h*(yold**2))/(1-2*h*yold)
        err = abs(ynew-yold)
        yold = ynew
    y.append(ynew)
    t.append(t[-1]+h)

T = time.clock() - T0   # Stop timing


print "Exact value y(%f) = %f " % (t[N], 1/(1-t[N]))
print "Value given by approx y_n(%f) = %f" % (t[N], y[N])
print "The error = y-y_n = %f " % (abs(1/(1-t[N]) - y[N]))
print "Time taken = %f " % (T)

yexact = [1.0/(1.0-x) for x in t]

plt.figure(); 
plt.plot(t,y,'r+',t,yexact,'b-.');
plt.xlabel('Time')
plt.ylabel('Solution')
plt.legend(('Backward Euler', 'Exact Solution'))
plt.title('Numerical solution of dy/dt=y^2')
plt.show()


