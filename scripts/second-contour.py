from numpy import *
from pylab import *

def wave(x,y,t=15.,a=36.,c=1.,theta=pi/4.):
    k0 = sqrt(2.)*pi/2.
    sint = sin(theta)
    cost = cos(theta)
    return ( exp(1.j*(k0*(x*cost-y*sint)-c*t)) * 
             exp(-((x*cost-y*sint)-c*t)**2/a) * 
             exp(-(x*sint+y*cost)**2/a) +
             exp(1.j*(k0*(x*cost+y*sint)-c*t)) * 
             exp(-((x*cost+y*sint)-c*t)**2/a) * 
             exp(-(x*sint-y*cost)**2/a) )

delta = 0.1
x = np.arange(-12.0, 12.0, delta)
y = np.arange(-12.0, 12.0, delta)
X, Y = np.meshgrid(x, y)

figure(figsize=(14,4))
suptitle('Contour plot of two interfering wave packets',
            fontsize=15)

ax = subplot(131)
#ax.text(200,220,"t=-15s")
text(0.9, 0.9,"t=-15s", 
        ha="center", va="center", transform=ax.transAxes)
xlabel("x")
ylabel("y")
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
Z = wave(X,Y,t=-15.)
contour(Z*conj(Z))

ax = subplot(132)
text(0.9, 0.9,"t=0s", 
        ha="center", va="center", transform=ax.transAxes)
xlabel("x")
ylabel("y")
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
Z = wave(X,Y,t=0)
contour(Z*conj(Z))

ax = subplot(133)
text(0.1, 0.9,"t=15s", 
        ha="center", va="center", transform=ax.transAxes)
xlabel("x")
ylabel("y")
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
Z = wave(X,Y,t=15.)
c = contour(Z*conj(Z))

cbar = colorbar(c)
cbar.set_ticks([])

legend()

show()
