#!/usr/bin/env python

# parse arguments
from switch import *
from sys import argv,stdin
expr = argv[1] if len(argv) > 1 else "p"
file = argv[2] if len(argv) > 2 else None

# import libraries
import matplotlib
matplotlib.use('macosx' if file is None else 'PDF')

import pyfits
import numpy as np
import scipy.stats
from pylab import *

# load trajectory data
hdu = pyfits.open('smpout.fit')

phi = hdu[0].data
chi = hdu[1].data
pip = hdu[2].data
pic = hdu[3].data
a   = hdu[4].data
p   = hdu[5].data

# load PDF data, J = P(t)/P(0)
pdf = pyfits.open('derivs.fit')

J_phi = pdf[0].data
J_chi = pdf[1].data
J_pip = pdf[2].data
J_pic = pdf[3].data
J_a   = pdf[4].data
J_p   = pdf[5].data

# variable being plotted
for case in switch(expr):
	if case("phi"): q = a*phi; break
	if case("chi"): q = a*chi; break
	if case("pip"): q = pip/a; break
	if case("pic"): q = pic/a; break
	if case("lna"): q = np.log(a); break
	if case("p"):   q = -p/6.0; break
	if case("lnF(phi)"): q = -np.log(1.0/np.square(J_phi) + 1.0e-16)/2.0; break
	if case("lnF(chi)"): q = -np.log(1.0/np.square(J_chi) + 1.0e-16)/2.0; break
	if case("lnF(pip)"): q = -np.log(1.0/np.square(J_pip) + 1.0e-16)/2.0; break
	if case("lnF(pic)"): q = -np.log(1.0/np.square(J_pic) + 1.0e-16)/2.0; break
	if case("lnF(a)"):   q = -np.log(1.0/np.square(J_a)   + 1.0e-16)/2.0; break
	if case("lnF(p)"):   q = -np.log(1.0/np.square(J_p)   + 1.0e-16)/2.0; break

np.putmask(q, np.isnan(q), 0.0)
np.putmask(q, np.isinf(q), 0.0)

ny,nx = q.shape

# data extent and meshes
x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = x0 + arange(nx)*dx
y = y0 + arange(ny)*dy

X, Y = meshgrid(x, y)

hdu.close(); pdf.close()

# high-pass filter in parameter direction?
for case in switch(expr):
	if case("never"):
		for i in range(0,ny):
			q[i,:] = q[i,:] - sum(q[i,:])/nx
		break
	if case("lna", "p"):
		w = kaiser(nx/4+1, 12.0)
		v = convolve(np.ones(nx), w, mode='same')
		for i in range(0,ny):
			q[i,:] = q[i,:] - convolve(q[i,:], w, mode='same')/v
		break

# smooth in time direction?
for case in switch(expr):
	if case("never"):
		w = kaiser(ny/4+1, 12.0)
		v = convolve(np.ones(ny), w, mode='same')
		for j in range(0,nx):
			q[:,j] = convolve(q[:,j], w, mode='same')/v
		break

# smooth in parameter direction?
for case in switch(expr):
	if case("never"):
		w = kaiser(17, 12.0)
		v = convolve(np.ones(nx), w, mode='same')
		for i in range(0,ny):
			q[i,:] = convolve(q[i,:], w, mode='same')
		break

# remap extreme values
b = max(-q.min(), q.max())
for case in switch(expr):
	if case("phi","chi","pip", "pic"): Q = q/b; break
	if case(): Q = q/sqrt(b*b/16.0 + q*q); break

# colormap for plotting
for case in switch(expr):
	if case("lnF(phi)", "lnF(chi)", "lnF(pip)", "lnF(pic)", "lnF(a)", "lnF(p)"):
		gradient = ["plum", "blue", "palegreen", "red", "white"]; break
	if case():
		gradient = ["blue", "white", "red"]; break

# create the figure
figure(figsize=(10,8), frameon=False)
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", gradient)
imshow(Q, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=dx/dy, cmap=c, vmin=-1.0, vmax=1.0, interpolation='none')
#colorbar(shrink=0.585)
#contour(X, Y, R, 32, cmap=cm.jet)
#contour(X, Y, q, [0.0], colors='white', linewidths=2.0)

data3d = np.load('data/slice-1875.npz')
plot(data3d['x'], 50.0+data3d['dlna'], 'k-')

xlim([x[0],x[-1]])
#xticks(arange(5)-11.0)

xlabel('$\\alpha \equiv \ln(\chi_0/\phi_0)/\mu_0T$')
ylabel('$\\tau$')

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=True)
show()
