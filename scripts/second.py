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

if __name__ == "__main__":
    inits = linspace(-10,-30,50)
    inits = append(inits, linspace(10.5,30.5,50))
    Second(save=True,plot=True, inits=inits,t0=-15.,tf=15.)
