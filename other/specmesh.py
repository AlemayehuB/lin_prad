#!/usr/bin/python

import math
import numpy as np
from scipy.fftpack import fftn, ifftn
from scipy.stats import norm

meshfile = 'specmesh'

Ncell = 128
headstr = "# Ncell = %d\n" % Ncell

Ntaper = 4
headstr = "# Ntaper = %d\n" % Ntaper

B0 = 2.0E+04
headstr += "# B0 = %12.5E\n" % B0

plidx = -1.5
headstr += "# plidx = %12.5E\n" % plidx

nlo = 4.0
nhi = 40.0
headstr += "# nlo = %12.5E ; nhi = %12.5E\n" % (nlo, nhi)

rng_seed = 32344544
headstr += "# rng_seed = %d\n" % rng_seed

eps = 1.0 / Ntaper
plnrm = (plidx + 1) * B0**2 / (4096.0 * math.pi**9 * ( nhi**(plidx+1) - nlo**(plidx+1) ))
dvol = (2* math.pi)**3
delta = 1.0 / Ncell
NDFT = Ncell
HNDFT = NDFT/2
np.random.seed(rng_seed)

ftarr = np.zeros((3, NDFT, NDFT, NDFT), dtype=np.complex_)
kn = np.zeros(3)
zr = np.zeros(3)

print "Making Fourier-space array..."
for n1 in range(HNDFT):
    print n1
    for n2 in range(-HNDFT, HNDFT):
        for n3 in range(-HNDFT, HNDFT):
            if n1 == 0 and n2 == 0 and n3 == 0: continue

            kn[:] = (n1, n2, n3)
            magkn = math.sqrt(kn.dot(kn))
            if magkn>=nlo and magkn<=nhi:
                std = math.sqrt(plnrm * magkn ** (plidx-4))
                rv = norm.rvs(loc=0.0, scale=std, size=6)
                Abuf = np.array([rv[0] + 1j*rv[1], rv[2] + 1j*rv[3], rv[4] + 1j*rv[5]])
                ftarr[:, n1%NDFT, n2%NDFT, n3%NDFT] = Abuf
                ftarr[:, -n1%NDFT, -n2%NDFT, -n3%NDFT] = Abuf.conj() # Yes, I know, when n1==0 this happens twice.  It's harmless.

print "...done.  Transforming..."
Aarr = fftn(ftarr, axes=(1,2,3)).real.astype(float) * dvol

print "...done.  Tapering..."
fbuf = Aarr.reshape(3, NDFT**3)
bases = np.array([[NDFT**2,NDFT,1], [NDFT,1,NDFT**2], [1,NDFT**2,NDFT]])
for idir in range(3):
    for i0 in range(Ntaper):
        for i1 in range(NDFT):
            for i2 in range(NDFT):
                idx = i0 * bases[idir,0] + i1 * bases[idir,1] + i2 * bases[idir,2]
                fbuf[:,idx] *= i0*eps
                idx = (NDFT-1-i0) * bases[idir,0] + i1 * bases[idir,1] + i2 * bases[idir,2]
                fbuf[:,idx] *= i0*eps

del(fbuf)
del(ftarr)
print "...done.  Taking curl of A..."

b = np.zeros(3)
j = np.zeros(3)
Barr = np.zeros((3, Ncell, Ncell, Ncell))

def bc(ix,iy,iz,ind, arr):
    if ix < 0 or ix > Ncell-1 or iy < 0 or iy > Ncell-1 or iz < 0 or iz > Ncell-1:
        return 0.0
    else:
        return arr[ind,ix,iy,iz]

for ix in range(Ncell):
    for iy in range(Ncell):
        for iz in range(Ncell):
            Barr[:,ix, iy, iz] = \
              ((bc(ix,iy+1,iz,2, Aarr)-bc(ix,iy-1,iz,2, Aarr)) - (bc(ix,iy,iz+1,1, Aarr)-bc(ix,iy,iz-1,1, Aarr)),\
               (bc(ix,iy,iz+1,0, Aarr)-bc(ix,iy,iz-1,0, Aarr)) - (bc(ix+1,iy,iz,2, Aarr)-bc(ix-1,iy,iz,2, Aarr)),\
               (bc(ix+1,iy,iz,1, Aarr)-bc(ix-1,iy,iz,1, Aarr)) - (bc(ix,iy+1,iz,0, Aarr)-bc(ix,iy-1,iz,0, Aarr)) )
            Barr[:,ix, iy, iz] /= 2 * delta

print "...done.  Output time."

mf = open(meshfile, 'w')
mf.write(headstr)
mf.write('# Cells per dimension = %d\n' % Ncell)
divB = 0.0
RMSB = 0.0
nc = 0
for ix in range(Ncell):
    for iy in range(Ncell):
        for iz in range(Ncell):

            b = Barr[:,ix, iy, iz]
            j[:] = ((bc(ix,iy+1,iz,2, Barr)-bc(ix,iy-1,iz,2, Barr)) - (bc(ix,iy,iz+1,1, Barr)-bc(ix,iy,iz-1,1, Barr)),\
                    (bc(ix,iy,iz+1,0, Barr)-bc(ix,iy,iz-1,0, Barr)) - (bc(ix+1,iy,iz,2, Barr)-bc(ix-1,iy,iz,2, Barr)),\
                    (bc(ix+1,iy,iz,1, Barr)-bc(ix-1,iy,iz,1, Barr)) - (bc(ix,iy+1,iz,0, Barr)-bc(ix,iy-1,iz,0, Barr)))
            j /= 2*delta

            db =   bc(ix+1,iy,iz,0, Barr)-bc(ix-1,iy,iz,0, Barr) + \
                   bc(ix,iy+1,iz,1, Barr)-bc(ix,iy-1,iz,1, Barr) + \
                   bc(ix,iy,iz+1,2, Barr)-bc(ix,iy,iz-1,2, Barr)
            nrm = max(abs(bc(ix+1,iy,iz,0, Barr)), abs(bc(ix-1,iy,iz,0, Barr)), \
                      abs(bc(ix,iy+1,iz,0, Barr)), abs(bc(ix,iy-1,iz,0, Barr)), \
                      abs(bc(ix,iy,iz+1,0, Barr)), abs(bc(ix,iy,iz-1,0, Barr)) )
            b2 = b.dot(b)
            if nrm > 0.0:
                divB += abs(db) / nrm
                RMSB += b2
                nc += 1

            buf = "%4d %4d %4d " % (ix, iy, iz)
            for bval in b: buf += "%17.10E " % (bval)
            for jval in j: buf += "%17.10E " % (jval)

            mf.write("%s\n" % buf)

mf.close()

RMSB = math.sqrt(RMSB / nc)
divB /= nc
print "RMS B = %12.5E G ; Normalized UDD Divergence = %12.5E G, # of cells = %d." % (RMSB, divB, nc)
