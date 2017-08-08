#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math
import sys
import re
from scipy.fftpack import fftn, ifftn
import rad_ut as ru
from rad_ut import *

# Margin
marg = 0.98

# Data file
dfile = sys.argv[1]

# Output file
ofile = sys.argv[2]


# Proton limiter
plimit = -1
if len(sys.argv) > 3:
    plimit = int(sys.argv[3])

# Resolution of plot
nbins = 128
if len(sys.argv) > 4:
    nbins = int(sys.argv[4])
nel = nbins**2
    
# GS iteration tolerance
tol_iter = 1.0E-4

# Maximum number of GS iterations
max_iter = 4000


################################################################################

out = open(ofile, 'w')
fd = open(dfile, 'r')

out.write("# Input filename: %s\n" % dfile)

line = fd.readline()
while not re.search('^# Columns:', line):
    
    if re.search('^# Tkin:', line):
        Tkin = float(line.split()[2])
    
    if re.search('^# rs:', line):
        rs = float(line.split()[2])
    
    if re.search('^# ri:', line):
        ri = float(line.split()[2])
    
    if re.search('^# raperture:', line):
        rap = float(line.split()[2])
    
    out.write(line)
    line = fd.readline()
    
while re.search('^#', line): line = fd.readline()

radius = rap * rs / ri  # radius of undeflected image of aperture at screen
ru.dmax = marg * radius / math.sqrt(2.0)
ru.delta = 2.0 * ru.dmax / nbins
delta_i = ru.delta * ri / rs

e = 4.8032E-10   # Statcoul
mp = 1.6726E-24  # g
c = 2.9979E+10   # cm/s
ergperMeV = 1.6022E-06 # 1 MeV / 1 erg
v = math.sqrt(2*(Tkin*ergperMeV)/mp)
Bconst = mp * c * v / (e * (rs-ri))

C = np.zeros((nbins,nbins))
J = np.zeros((nbins,nbins))
JS = np.zeros((nbins,nbins))
Lam = np.zeros((nbins,nbins))
LamS = np.zeros((nbins,nbins))
ExpLam = np.zeros((nbins,nbins))
Src = np.zeros((nbins,nbins))
Bperp = np.zeros((nbins,nbins,2))
BperpR = np.zeros((nbins,nbins,2))
BperpS = np.zeros((nbins,nbins,2))
deltaX = np.zeros((nbins,nbins,2))
deltaXR = np.zeros((nbins,nbins,2))
deltaXS = np.zeros((nbins,nbins,2))

buf = "# nbins = %d ; delta = %12.5E\n" % (nbins, ru.delta)
out.write(buf)

print "Reading %s, binning %d x %d..." % (dfile, nbins, nbins)

nprot = 0
while line:
    nprot += 1
    nx = float(line.split()[0])
    ny = float(line.split()[1])
    xx= float(line.split()[3])
    yy= float(line.split()[4])
    jj = float(line.split()[8])
    b0 = float(line.split()[9])
    b1 = float(line.split()[10])
    idx = vec2idx((xx,yy))
    i = idx[0] ; j = idx[1]
    
    if (xx + ru.dmax)/ru.delta >= 0 and i < nbins and (yy + ru.dmax)/ru.delta >= 0 and j < nbins:
        C[i,j] += 1
        Bperp[i,j,0] += b0
        Bperp[i,j,1] += b1
        deltaX[i,j,0] += (xx - nx*rs)
        deltaX[i,j,1] += (yy - ny*rs)
        J[i,j] += jj

    if nprot == plimit: break
    line = fd.readline()
    
fd.close()

buf = "Read %d protons, of which %d in square window." % (nprot, C.sum())
print "done. " + buf
out.write("# %s\n" % buf)

print "Min, max, mean pixel counts, and delta:"
print C.min(), C.max(), C.mean(), ru.delta

avg_fluence = nprot / (math.pi * radius**2)

for i in range(nbins):
    for j in range(nbins):
        try: 
            Bperp[i,j,:] /= C[i,j]
            deltaX[i,j,:] /= C[i,j]
            J[i,j] /= C[i,j]
            Lam[i,j] = 2.0 * ( 1.0 - math.sqrt(avg_fluence * ru.delta**2/C[i,j]) )
            ExpLam[i,j] = math.exp(Lam[i,j])
            Src[i,j] = Lam[i,j] * ExpLam[i,j]
        except ZeroDivisionError:
            print "Zero pixel, will screw everything up."
            raise ValueError
        
print "Gauss-Seidel Iteration..."

def D(i,j):
    d = -2.0 * ExpLam[i,j] - 0.5 * ( bc_enforce_N(ExpLam, i+1,j) + \
                                     bc_enforce_N(ExpLam, i-1,j) + \
                                     bc_enforce_N(ExpLam, i,j+1) + \
                                     bc_enforce_N(ExpLam, i,j-1) )
    return d

def O(i,j, x):
    a = 0.5 * ( bc_enforce_D(x, i+1,j) * (bc_enforce_N(ExpLam, i+1,j) + ExpLam[i,j]) +\
                bc_enforce_D(x, i-1,j) * (bc_enforce_N(ExpLam, i-1,j) + ExpLam[i,j]) +\
                bc_enforce_D(x, i,j+1) * (bc_enforce_N(ExpLam, i,j+1) + ExpLam[i,j]) +\
                bc_enforce_D(x, i,j-1) * (bc_enforce_N(ExpLam, i,j-1) + ExpLam[i,j]) )
    return a

# Initial guess by FFT Poisson solve
phi = solve_poisson(Lam)

# Iterate to solution
GS = Gauss_Seidel(phi, D, O, Src, talk=20, tol=tol_iter, maxiter=max_iter)

phi *= ru.delta**2

for i in range(nbins):
    for j in range(nbins):
        deltaXR[i,j] = -gradient(phi, (i,j))
        BperpR[i,j,0] = Bconst * deltaXR[i,j,1]
        BperpR[i,j,1] = -Bconst * deltaXR[i,j,0]
        x = idx2vec((i,j))
        x = x + deltaXR[i,j]
        idx = vec2idx(x)
        BperpS[i,j,:] = Bperp[idx[0]%nbins, idx[1]%nbins, :]
        deltaXS[i,j,:] = deltaX[idx[0]%nbins, idx[1]%nbins, :]
        JS[i,j] = J[idx[0]%nbins, idx[1]%nbins]
        LamS[i,j] = Lam[idx[0]%nbins, idx[1]%nbins]

EB = BperpS.flatten().dot(BperpS.flatten()) * delta_i**2
EBR = BperpR.flatten().dot(BperpR.flatten()) * delta_i**2
L2B = fnorm(BperpR-BperpS)/fnorm(BperpS)

buf = "EB = %12.5E ; EB_Reconstructed = %12.5E" % (EB, EBR)
print "...done.  " + buf
out.write("# %s\n" % buf)

buf = "Relative L2 norm of (reconstructed B - actual B) = %12.5E" % L2B
print buf
out.write("# %s\n" % buf)

print "Output time!"

out.write("#L2 norm of residual = %12.5E ;  Number of Gauss-Seidel iterations = %d\n" % GS)

out.write('# Column 1: x index\n')
out.write('# Column 2: y index\n')
out.write('# Column 3: x\n')
out.write('# Column 4: y\n')
out.write('# Column 5: Bx (actual)\n')
out.write('# Column 6: By (actual)\n')
out.write('# Column 7: Bx (solver)\n')
out.write('# Column 8: By (solver)\n')
out.write('# Column 9: Delta_x (actual)\n')
out.write('# Column 10: Delta_y (actual)\n')
out.write('# Column 11: Delta_x (solver)\n')
out.write('# Column 12: Delta_y (solver)\n')
out.write('# Column 13: Position-corrected Lambda\n')
out.write('# Column 14: Position-corrected current integral\n')

for i in range(nbins):
    for j in range(nbins):

        x = idx2vec((i,j))
        line = ("%3d  %3d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  " +\
               "%12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n") %\
               (i, j, x[0], x[1], \
                BperpS[i,j,0], BperpS[i,j,1], \
                BperpR[i,j,0], BperpR[i,j,1], \
                deltaXS[i,j,0], deltaXS[i,j,1], \
                deltaXR[i,j,0], deltaXR[i,j,1], \
                LamS[i,j], JS[i,j] )
                
        out.write(line)
out.close()



