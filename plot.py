#!/usr/bin/env python

import matplotlib
matplotlib.use('macosx')

import pyfits
import numpy as np
import scipy.stats
from pylab import *

hdu = pyfits.open('smpout.fit')

phi = hdu[0].data
chi = hdu[1].data
pip = hdu[2].data
pic = hdu[3].data
a   = hdu[4].data
p   = hdu[5].data

q = -p/6.0

np.putmask(q, np.isnan(q), 0.0)
np.putmask(q, np.isinf(q), 0.0)

ny,nx = q.shape

x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = (x0+dx) + arange(nx)*dx
y = (y0+dy) + arange(ny)*dy

X, Y = meshgrid(x, y)

hdu.close()

## destripe in time direction
for i in range(0,ny):
	q[i,:] = q[i,:] - sum(q[i,:])/nx

# smooth in time direction
w = kaiser(ny/4+1, 12.0)
#for j in range(0,nx):
#	q[:,j] = convolve(q[:,j], w, mode='same')

# smooth in parameter direction
w = kaiser(17, 12.0)
#for i in range(0,ny):
#	q[i,:] = convolve(q[i,:], w, mode='same')

# colormap for plotting
b = max(-q.min(), q.max())/4.0; Q = q #/sqrt(b*b + q*q)

# create the figure
figure(figsize=(10,8), frameon=False)

c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "red"])
#c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["plum", "blue", "palegreen", "red", "white"])
imshow(Q, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=dx/dy, cmap=c, vmin=-1.0, vmax=1.0, interpolation='none')
#colorbar(shrink=1.0)
#contour(X, Y, R, 32, cmap=cm.jet); contour(X, Y, q, [0.0], colors='white', linewidths=2.0)

data3d = np.load('data/slice-1875.npz')
plot(data3d['x'], 50.0+data3d['dlna'], 'k-')

xticks(arange(5)-11.0)
#yticks(arange(7)-3.0)

xlabel('$\\alpha \equiv \ln(\chi_0/\phi_0)/\mu_0T$')
ylabel('$\\tau$')

show()

#savefig("smpout.pdf", format="pdf")