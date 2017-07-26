#!/usr/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import poisson
import math
import sys
import re

# Data file
dfile = sys.argv[1]

fd = open(dfile, 'r')
line = fd.readline()
while re.search('^#', line):

    if re.search('^# nbins', line):
        nbins = int(line.split()[3])
        delta = float(line.split()[7])
        
    line = fd.readline()

X = np.zeros((nbins,nbins))
Y = np.zeros((nbins,nbins))
Br = np.zeros((nbins,nbins,2))
Bi = np.zeros((nbins,nbins,2))
BiMag = np.zeros((nbins,nbins))
BrMag = np.zeros((nbins,nbins))
dxr = np.zeros((nbins,nbins,2))
dxi = np.zeros((nbins,nbins,2))
dxiMag = np.zeros((nbins,nbins))
LamS = np.zeros((nbins,nbins))
JS = np.zeros((nbins,nbins))


while line:
    i = int(line.split()[0])
    j = int(line.split()[1])
    
    X[i,j] = float(line.split()[2])
    Y[i,j] = float(line.split()[3])
    
    Bi[i,j,0] = float(line.split()[4])
    Bi[i,j,1] = float(line.split()[5])
    BiMag[i,j] = 0.5*math.log10(Bi[i,j,0]**2 + Bi[i,j,1]**2)
    
    Br[i,j,0] = float(line.split()[6])
    Br[i,j,1] = float(line.split()[7])
    BrMag[i,j] = 0.5*math.log10(Bi[i,j,0]**2 + Bi[i,j,1]**2)

    
    line = fd.readline()

stretch = 13.1/10.2
plt.rc('text', usetex=True)

#################################################
fig  = plt.figure()
fig.set_figwidth(6.0 * stretch)
fig.set_figheight(6.0)

ax = fig.add_subplot(1,1,1)
strm = ax.streamplot(X[:,0], Y[0,:], Bi[:,:,0].T, Bi[:,:,1].T, color=BiMag.T, \
                      linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0)
fig.colorbar(strm.lines)
ax.set_title(r"Log True $B_\perp$ Projection (G cm)", fontsize=18)
ax.set_xlabel(r"X (cm)", fontsize=18)
ax.set_ylabel(r"Y (cm)", fontsize=18)
ax.tick_params(labelsize='large')

fig.savefig("B_True.png", format='png')

#################################################
fig.clf()
ax = fig.add_subplot(1,1,1)
strm = ax.streamplot(X[:,0], Y[0,:], Br[:,:,0].T, Br[:,:,1].T, color=BrMag.T, \
                      linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0)
fig.colorbar(strm.lines)
ax.set_title(r"Log Reconstructed $B_\perp$ Projection (G cm)", fontsize=18)
ax.set_xlabel(r"X (cm)", fontsize=18)
ax.set_ylabel(r"Y (cm)", fontsize=18)
ax.tick_params(labelsize='large')

fig.savefig("B_Reconstructed.png", format='png')

