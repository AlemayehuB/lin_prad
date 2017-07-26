#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.stats import poisson
import math
import sys
import re
import rad_ut as ru
from rad_ut import *

# Floor counts per bin
Cmin = 10.0

# Margin
marg = 0.98

# Proton limiter
plimit = 10000000

# Data file
dfile = sys.argv[1]
ofile = dfile + "_binned"

# Resolution of plot
nbins = 128
if len(sys.argv) > 2:
    nbins = int(sys.argv[2])

fd = open(dfile, 'r')
out = open(ofile, 'w')
line = fd.readline()

while not re.search('^# Tkin:', line):
    out.write(line)
    line = fd.readline()
Tkin = float(line.split()[2])

while not re.search('^# rs:', line):
    out.write(line)
    line = fd.readline()
rs = float(line.split()[2])

while not re.search('^# ri:', line):
    out.write(line)
    line = fd.readline()
ri = float(line.split()[2])

while not re.search('^# raperture:', line):
    out.write(line)
    line = fd.readline()
rap = float(line.split()[2])

while re.search('^#', line):
    out.write(line)
    line = fd.readline()

radius = rap * rs / ri  # radius of undeflected image of aperture at screen
dmax = marg * radius / math.sqrt(2.0)
delta = 2.0 * dmax / nbins
ru.dmax = dmax
ru.delta = delta


buf = "# nbins = %d ; delta = %12.5E\n" % (nbins, ru.delta)
out.write(buf)

e = 4.8032E-10   # Statcoul
mp = 1.6726E-24  # g
c = 2.9979E+10   # cm/s
ergperMeV = 1.6022E-06 # 1 MeV / 1 erg
v = math.sqrt(2*(Tkin*ergperMeV)/mp)
 
Bconst = delta * mp * c * v * rs / (math.pi * e * ri * (rs-ri) )


x = np.zeros((nbins+1, nbins+1))
y = np.zeros((nbins+1, nbins+1))
for i in range(nbins+1):
    xx = -dmax + i*delta
    for j in range(nbins+1):
        yy = -dmax + j*delta
        x[i,j] = xx
        y[i,j] = yy

C = np.zeros((nbins,nbins))
J = np.zeros((nbins,nbins))
Theta = np.zeros((nbins,nbins))
Diff = np.zeros((nbins,nbins))
DiffCts = np.zeros((nbins,nbins))
Fluct = np.zeros((nbins,nbins))
eta = np.zeros((nbins,nbins))
Bperp = np.zeros((nbins,nbins,2))
nprot = 0

xpt=[]
ypt=[]
while line:
    nprot += 1
    n1 = np.array([float(s) for s in line.split()[0:3]])
    n2 = np.array([float(s) for s in line.split()[5:8]])
    theta = math.acos(min(1.0, n1.dot(n2)/math.sqrt(n1.dot(n1)*n2.dot(n2))))
    xx= float(line.split()[3])
    yy= float(line.split()[4])
    jj = float(line.split()[8])
    i = int((xx + dmax)/delta)
    j = int((yy + dmax)/delta)
    if (xx + dmax)/delta >= 0 and i < nbins and (yy + dmax)/delta >= 0 and j < nbins:
        C[i,j] += 1
        J[i,j] += jj
#        eta[i,j] += 1.0 / (1.0 - 0.5 * jj)**2
        Theta[i,j] = theta
        xpt.append(xx)
        ypt.append(yy)

    if nprot == plimit: break
    line = fd.readline()

avg_fluence = nprot / (math.pi * radius**2)
im_fluence = C.sum() / (4 * dmax**2)
for i in range(nbins):
    for j in range(nbins):
        if C[i,j] != 0:
            J[i,j] /= C[i,j]  # average current in cell
#            eta[i,j] = avg_fluence * delta**2 * eta[i,j] / C[i,j] # average of predicted fluence in cell
            eta[i,j] = avg_fluence * delta**2 / ( 1.0 - 0.5*J[i,j] )**2

        if C[i,j] >= Cmin:
            Fluct[i,j] = 2.0 * ( 1.0 - math.sqrt(avg_fluence * delta**2/C[i,j]) )
            err = math.sqrt(avg_fluence * delta**2 ) / C[i,j]
            Diff[i,j] = ( Fluct[i,j] - J[i,j] ) / err
            DiffCts[i,j] = ( C[i,j] - eta[i,j] ) / math.sqrt(C[i,j])


Dpos = Diff[ C>=Cmin ]
Jpos = J[ C>=Cmin ]
Flpos = Fluct[ C>=Cmin ]
DiffCtspos = DiffCts[ C>=Cmin ]
etapos = eta[ C>=Cmin ]

print "Total protons in sample: %d" % nprot
print "Total sample Average Fluence: %12.5E ; Image Average Fluence: %12.5E" % (avg_fluence, im_fluence)

print "Mean J: %12.5E ; Std. Dev. J: %12.5E" % (Jpos.mean(), Jpos.std())
print "Max J: %12.5E ; Min J: %12.5E" % (Jpos.max(), Jpos.min())

print "Mean Fluct.: %12.5E ; Std. Dev. Fluct.: %12.5E" % (Flpos.mean(), Flpos.std())
print "Max Fluct: %12.5E ; Min Fluct: %12.5E" % (Flpos.max(), Flpos.min())

print "Mean difference: %12.5E ; Std. Dev. Difference: %12.5E" % (Dpos.mean(), Dpos.std())
print "Max difference: %12.5E ; Min difference: %12.5E" % (Dpos.max(), Dpos.min())

print "Mean  predicted counts per bin: %12.5E ; Std. Dev. predicted counts per bin: %12.5E" % (etapos.mean(), etapos.std())
print "Max predicted counts per bin: %d ; Min predicted counts per bin: %d" % (etapos.max(), etapos.min())

print "Mean counts per bin: %12.5E ; Std. Dev. Counts per bin: %12.5E" % (C.mean(), C.std())
print "Max counts per bin: %d ; Min counts per bin: %d" % (C.max(), C.min())
print "Number of bins with zero protons: %d" % (C.size - C[ C>0 ].size)
print "Number of bins with %d or fewer protons: %d" % (Cmin, C.size - C[ C>Cmin ].size)

print "Mean (Count-Pred)/err: %12.5E ; Std. Dev.: %12.5E" % (DiffCtspos.mean(), DiffCtspos.std())
print "Max (Count-Pred)/err: %12.5E ; Min (Count-Pred)/err: %12.5E" % (DiffCtspos.max(), DiffCtspos.min())


print "Total protons in image: %d ; Total predicted: %d" % (C.sum(), eta.sum())


out.write('# Column 1: x index\n')
out.write('# Column 2: y index\n')
out.write('# Column 3: x\n')
out.write('# Column 4: y\n')
out.write('# Column 5: Fluence Contrast\n')
out.write('# Column 6: Current Projection\n')
out.write('# Column 7: Fluence Contrast Difference/Noise\n')
out.write('# Column 8: Counts/Bin\n')
out.write('# Column 9: Predicted Counts/Bin\n')
out.write('# Column 10: Counts Difference/Noise\n')


for i in range(nbins):
    for j in range(nbins):

        xx = idx2vec((i,j))
        line = ("%3d  %3d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  " +\
               "%12.5E  %12.5E\n") %\
               (i, j, xx[0], xx[1], \
                Fluct[i,j], J[i,j], Diff[i,j], C[i,j], eta[i,j], DiffCts[i,j] )

        out.write(line)
out.close()



stretch = 13.0 / 9.0

fig  = plt.figure()
fig.set_figwidth(18.0 * stretch)
fig.set_figheight(12.0)

ax = fig.add_subplot(2,3,1)
p = ax.pcolormesh(x,y, Fluct, cmap=cm.afmhot, vmin=min(Fluct.min(), J.min()), vmax=max(Fluct.max(), J.max()))
ax.set_ylabel("Y (cm)")
plt.colorbar(p)
ax.set_title("Fluence Contrast")

ax = fig.add_subplot(2,3,2)
p = ax.pcolormesh(x,y, J, cmap=cm.afmhot, vmin=min(Fluct.min(), J.min()), vmax=max(Fluct.max(), J.max()))
plt.colorbar(p)
ax.set_title("Current Projection Function")

vm = max(abs(Diff.min()), abs(Diff.max()))
ax = fig.add_subplot(2,3,3)
p = ax.pcolormesh(x,y, Diff, cmap=cm.RdYlGn, vmax=vm, vmin=-vm)
plt.colorbar(p)
ax.set_title("Fluence Contrast Difference/Noise")


ax = fig.add_subplot(2,3,4)
p = ax.pcolormesh(x,y, C, cmap=cm.afmhot, vmin=min(eta.min(), C.min()), vmax=max(eta.max(), C.max()))
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
plt.colorbar(p)
ax.set_title("Counts/Bin")

ax = fig.add_subplot(2,3,5)
p = ax.pcolormesh(x,y, eta, cmap=cm.afmhot, vmin=min(eta.min(), C.min()), vmax=max(eta.max(), C.max()))
ax.set_xlabel("X (cm)")
plt.colorbar(p)
ax.set_title("Predicted Counts/Bin")

vm = max(abs(DiffCts.min()), abs(DiffCts.max()))
ax = fig.add_subplot(2,3,6)
p = ax.pcolormesh(x,y, DiffCts, cmap=cm.RdYlGn, vmax=vm, vmin=-vm)
ax.set_xlabel("X (cm)")
plt.colorbar(p)
ax.set_title("Counts Difference/Noise")


fig.savefig("radiography.png", format='png')

fig.clf()

fig.set_figwidth(6.0)
fig.set_figheight(6.0)
ax = fig.add_subplot(1,1,1)

ax.plot(xpt, ypt, '.', markersize = 0.01 * math.sqrt(1.0E+06/nprot))
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
ax.set_title("Proton Locations")
fig.savefig("radiography_points.png", format='png')



fig.clf()

fig.set_figwidth(6.0/0.85)
fig.set_figheight(6.0)
ax = fig.add_subplot(1,1,1)
p = ax.pcolormesh(x,y, Theta, cmap=cm.RdYlGn)
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
ax.set_title("Angular Deviation")
plt.colorbar(p)

fig.savefig("radiography_angdev.png", format='png')
