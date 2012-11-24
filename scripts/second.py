import sympy as sp

from numpy import *
from bohm import Bohm

class Second(Bohm):
    def psi(self):
        return self.wavefunction(self.y,self.t)

    def wavefunction(self,y,t):
        # initial conditions
        a = 3.
        c = 1.
        k0 = sqrt(2.)*pi/2.
        x = 1.
        sint = sin(pi/4.)
        cost = cos(pi/4.)
        return ( sp.exp(1.j*(k0*(x*cost-y*sint)-c*t)) * 
                 sp.exp(-((x*cost-y*sint)-c*t)**2/4/a**2) * 
                 sp.exp(-(x*sint+y*cost)**2/4/a**2) +
                 sp.exp(1.j*(k0*(x*cost+y*sint)-c*t)) * 
                 sp.exp(-((x*cost+y*sint)-c*t)**2/4/a**2) * 
                 sp.exp(-(x*sint-y*cost)**2/4/a**2) )

if __name__ == '__main__':
    inits = array([x for x in random.uniform(-30,30,20)])
    Second(plot=True, inits=inits,t0=-5.,tf=5.,dt=1/10.)
