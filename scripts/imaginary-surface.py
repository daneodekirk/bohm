from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from pylab import *

def psi_y(y,t=0.):
    Y  = 1.
    s0 = .2
    k  = .1
    x  = 0.
    st  = s0*(1+ (1j*t)/(2*s0**2))
    N  = (2*pi*st**2)**-.25
    return N*exp((-((y-Y)**2)/(4*s0*st)) + 1j*(k*x-((t*k**2)/2)))

t = arange(-0,1.1,.01)
y = arange(-13,13,.01)
Y,T = meshgrid(y, t)
Z = psi_y(Y,T)+psi_y(-Y,T)

fig = figure(figsize=(17,10))
suptitle("Imaginary surface over time of the wave function for a double slit system",
            fontsize=20)

ax = fig.gca(projection="3d")
ax.set_zlim([-10,10])

surf = ax.plot_surface(Y, T, imag(Z), cmap=cm.jet, linewidth=0, antialiased=False)
ax.set_xlabel("y",fontsize=20)
ax.set_ylabel("t",fontsize=20)
ax.set_zlabel(r"$\psi$",fontsize=30)
ax.view_init(elev=36, azim=63) 

cbar = colorbar(surf)
cbar.set_label("Surface gradient", fontsize=20)

show()
