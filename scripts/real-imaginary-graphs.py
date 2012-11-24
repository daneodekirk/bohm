from scipy import *
from numpy import *
from pylab import *

def psi_y(y,t=1.):
    Y  = 1.
    s0 = .2
    k  = .1
    x  = 0.
    st  = s0*(1+ (1j*t)/(2*s0**2))
    N  = (2*pi*st**2)**-.25
    return N*exp((-((y-Y)**2)/(4*s0*st)) + 1j*(k*x-((t*k**2)/2)))

y = arange(-13,13,.01)
Z = psi_y(y)+psi_y(-y)

figure(figsize=(12,10))

suptitle("Real and imaginary parts of the wave function for the double slit system", 
            fontsize=20)

subplot(211)
xlabel("y",fontsize=20)
ylabel(r"$\psi$",fontsize=30)
plot(y,real(Z), "b", label="Real part")
legend()

subplot(212)
xlabel("y",fontsize=20)
ylabel(r"$\psi$",fontsize=30)
plot(y,imag(Z), "r", label="Imaginary part")
legend()

show()
