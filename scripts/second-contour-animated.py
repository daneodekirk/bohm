import os
import time

from numpy import *
from pylab import *

cwd = os.getcwd() + "/temp"

def wave(x,y,t=15.):
    k0 = sqrt(2.)*pi/2.
    a = 3.
    c = 1.
    sint = sin(pi/4.)
    cost = cos(pi/4.)
    return ( exp(1.j*(k0*(x*cost-y*sint)-c*t))*5.*sqrt(pi) * 
             exp(-((x*cost-y*sint)-c*t)**2./4./a**2.)/a*5.*sqrt(pi) * 
             exp(-(x*sint+y*cost)**2./4./a**2.)/a +
             exp(1.j*(k0*(x*cost+y*sint)-c*t))*5.*sqrt(pi) * 
             exp(-((x*cost+y*sint)-c*t)**2./4./a**2.)/a*5.*sqrt(pi) * 
             exp(-(x*sint-y*cost)**2./4./a**2.)/a )

delta = 0.1
x = np.arange(-12.0, 12.0, delta)
y = np.arange(-12.0, 12.0, delta)
X, Y = np.meshgrid(x, y)

t = -25.
while t<25.:
    clf()
    Z = wave(X,Y,t)
    contour(Z*conj(Z))
    t+=1.
    pause(.001)
    if t == 25:
        t = -25.
