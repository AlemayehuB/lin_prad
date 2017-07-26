#!/usr/bin/python

import math
import numpy as np
from scipy.fftpack import fftn, ifftn
from scipy.stats import norm
from scipy.linalg import cho_solve, cho_factor

meshfile = 'specmesh_GP'

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

SE_sigma = 2.0
headstr += "# SE_sigma = %12.5E\n" % SE_sigma

GP_stencil = np.array([[0,0,0], [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]], dtype=np.int_)
for p in GP_stencil: 
    headstr += "# GP Stencil point: "
    for i in p: headstr += "%d " % i
    headstr += "\n"

eps = 1.0 / Ntaper
plnrm = (plidx + 1) * B0**2 / (4096.0 * math.pi**9 * ( nhi**(plidx+1) - nlo**(plidx+1) ))
dvol = (2* math.pi)**3 
delta = 1.0 / Ncell
NDFT = Ncell
HNDFT = NDFT/2
np.random.seed(rng_seed)
sm2 = 1.0/SE_sigma**2
sm4 = sm2**2
sm6 = sm2**3

Nstencil = len(GP_stencil)
Covariance = np.zeros((Nstencil, 3, Nstencil, 3))
for p in range(Nstencil):
    for pp in range(Nstencil):
        dx = (GP_stencil[p] - GP_stencil[pp]).astype(float)
        dx2 = dx.dot(dx)
        e = math.exp(-0.5*dx2 / sm2)
        for i in range(3):
            Covariance[p,i,pp,i] = 2*sm2 - dx2*sm4
            for ii in range(3):
                Covariance[p,i,pp,ii] += sm4 * dx[i] * dx[ii]
        Covariance[p,:,pp,:] *= e
        
CovM = Covariance.reshape(3*Nstencil, 3*Nstencil)
CCF = cho_factor(CovM)

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
        
def bcvec(ix,iy,iz, arr):
    if ix < 0 or ix > Ncell-1 or iy < 0 or iy > Ncell-1 or iz < 0 or iz > Ncell-1:
        return np.zeros(3)
    else:
        return arr[:,ix,iy,iz]


for ix in range(Ncell):
    if ix%8 == 0: print ix
    for iy in range(Ncell):
        for iz in range(Ncell):
            Barr[:,ix, iy, iz] = \
              ((bc(ix,iy+1,iz,2, Aarr)-bc(ix,iy-1,iz,2, Aarr)) - (bc(ix,iy,iz+1,1, Aarr)-bc(ix,iy,iz-1,1, Aarr)),\
               (bc(ix,iy,iz+1,0, Aarr)-bc(ix,iy,iz-1,0, Aarr)) - (bc(ix+1,iy,iz,2, Aarr)-bc(ix-1,iy,iz,2, Aarr)),\
               (bc(ix+1,iy,iz,1, Aarr)-bc(ix-1,iy,iz,1, Aarr)) - (bc(ix,iy+1,iz,0, Aarr)-bc(ix,iy-1,iz,0, Aarr)) )
            Barr[:,ix, iy, iz] /= 2 * delta

print "...done.  Computing weight vectors..."
Warr = np.zeros((Ncell, Ncell, Ncell, 3*Nstencil))
RHS = np.zeros((Nstencil, 3))
for ix in range(Ncell):
    if ix%8 == 0: print ix
    for iy in range(Ncell):
        for iz in range(Ncell):
            for ist in range(Nstencil):
                t = np.array([ix,iy,iz]) + GP_stencil[ist]
                RHS[ist,:] = bcvec(t[0], t[1], t[2], Barr)
            Warr[ix,iy,iz,:] = cho_solve(CCF, RHS.flatten())

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
            for bval in b: buf += "%22.16E " % (bval)
            for jval in j: buf += "%22.16E " % (jval)
            for wval in Warr[ix,iy,iz,:]: buf += "%22.16E " % (wval)

            mf.write("%s\n" % buf)

mf.close()

RMSB = math.sqrt(RMSB / nc)
divB /= nc
print "RMS B = %12.5E G ; Normalized UDD Divergence = %12.5E G, # of cells = %d." % (RMSB, divB, nc)

