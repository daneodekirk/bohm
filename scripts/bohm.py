import time
import os
import numpy as np
import sympy as sp

from pylab import *
from random import getrandbits

class Bohm:
    def __init__(self,particles=200,plot=False,save=False,**kwargs):
        self.start = time.time()

        self.slitdistance = kwargs.get('slitdistance',1.)
        self.slitwidth    = kwargs.get('slitwidth',.2)
        self.velx         = kwargs.get('velx',.1)
        self.dt           = kwargs.get('dt',1/100.)
        self.save         = save

        self.y = sp.var("y")
        self.t = sp.var("t")
        self.x  = y*t

        self.st = self.slitwidth*(1+ (1j*t)/(2*self.slitwidth**2))
        self.N  = (2*pi*self.st**2)**-.25

        inits = np.array([x for x in np.random.uniform(-1.5, 1.5, particles*2) 
                                    if x < -.5 or x > .5][:particles])

        self.inits     = kwargs.get('inits', inits)
        self.t0        = kwargs.get('t0',0)
        self.tf        = kwargs.get('tf',2)

        self.PSI = self.psi()
        # turn the symoblic functions into callable python functions
        self.wave         = sp.lambdify([y,t],self.PSI,"numpy")
        self.del_wave     = sp.lambdify([y,t],self.PSI.diff(y),"numpy")

        if plot: self.plot()

    def rk4(self,x, h, y, f):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5*h, y + 0.5*k1)
        k3 = h * f(x + 0.5*h, y + 0.5*k2)
        k4 = h * f(x + h, y + k3)
        return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

    def psi(self):
        return (self.wavefunction(self.x,self.y,self.t) 
                    + self.wavefunction(self.x,-self.y,self.t))

    def wavefunction(self,x,y,t):
        return self.N*sp.exp((-((y-self.slitdistance)**2.)/(4.*self.slitwidth*self.st)) 
                                + 1.j*(self.velx*x-((t*self.velx**2.)/2.)))

    def trajectories(self,t,state):
        y,vel = state
        return imag(self.del_wave(y,t)/self.wave(y,t))

    def plot(self):
        figure(figsize=(12,10))
        suptitle("Bohmian trajectories for %d particles" % len(self.inits), fontsize=20)
        xlabel("t")
        ylabel("y")
        for index,y0 in enumerate(self.inits):
            t = self.t0
            state = np.array([y0,0])
            while t < self.tf:
                t, state = self.rk4(t,self.dt,state,self.trajectories)
                plot(t,state[0],"o",markersize=1)

        print "Completed in %s seconds" % str(time.time() - self.start)
        if self.save: 
            filename = ("%s/Desktop/%d-bohmian-trajectories-%s.png" % 
                            (os.path.expanduser("~"), len(self.inits), os.urandom(8).encode('hex')))
            savefig(filename)
            print "Saved to %s" % filename
        show()

if __name__ == '__main__':
    Bohm(plot=True,save=True)
