#!/usr/bin/python
'''
 The plots  'fluence contrast', 'current projection function',
'Fluence Contrast Difference/Noise','Counts/Bin','Predicted Counts/Bin',
'Counts Difference/Noise',  'Angular Deviation',and 'Proton Locations' are all
created from the input file of reconstruct.py.
'''
from re import match
from sys import argv
import math

import rad_ut as ru

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# Data file
dfile = argv[1]
data = open(dfile, 'r')

# Output file
ofile = dfile + "_binned"
out = open(ofile, 'w')

# Resolution of plot
num_bins = 128
if len(argv) > 2: num_bins = int(argv[2])

#############Parsing Data File###################################################
line=data.readline()
while not match('# Column', line):

    if match('# Tkin:', line):
        Tkin = float(line.split()[2])  # Kinetic energy

    elif match('# rs:', line):
        rs = float(line.split()[2])  # Length from implosion to screen

    elif match('# ri:', line):
        ri = float(line.split()[2])  # Length from implosion to interaction region

    elif match('# raperture:', line):
        rap = float(line.split()[2])  # Aperture of the cone
    out.write(line)
    line=data.readline()

# The loop reads each line of the input until the data is reached.
while match('#', line): line = data.readline()

#############Setting Varibales##################################################
# Floor counts per bin
Cmin = 10.0

# Margin
marg = 0.98

# Proton limiter
plimit = 10000000

radius = rap * rs / ri  # radius of the undeflected image on the screen
dmax = marg * radius / math.sqrt(2.0)
delta = 2.0 * dmax / num_bins

# Recovered protons that hit the collimated square image
rec_prot=0

x = np.zeros((num_bins+1, num_bins+1))
y = np.zeros((num_bins+1, num_bins+1))
for i in range(num_bins+1):
    xx = -dmax + i*delta
    for j in range(num_bins+1):
        yy = -dmax + j*delta
        x[i,j] = xx
        y[i,j] = yy
# Pixel Count
count = np.zeros((num_bins , num_bins))
# Current integrals
cur_int= np.zeros((num_bins, num_bins))
# Angular Deviation
Theta = np.zeros((num_bins,num_bins))
# Difference
Diff = np.zeros((num_bins,num_bins))
# Counts Difference/Noise
DiffCts = np.zeros((num_bins,num_bins))
# Fluence Contrast
Fluct = np.zeros((num_bins,num_bins))
# Predicted Counts per Bin
eta = np.zeros((num_bins,num_bins))
# B Field Integral
B= np.zeros((num_bins,num_bins,2))
# Final X location at screen (cm)
xpt=[]
# Final Y location at screen (cm)
ypt=[]
################################################################################
while line:
    rec_prot += 1
    n1 = np.array([float(s) for s in line.split()[0:3]])
    n2 = np.array([float(s) for s in line.split()[5:8]])
    theta = math.acos(min(1.0, n1.dot(n2)/math.sqrt(n1.dot(n1)*n2.dot(n2))))
    xx= float(line.split()[3])  # Final X location at screen (cm)
    yy= float(line.split()[4])  # Final Y location at screen (cm)
    current= float(line.split()[8]) # Current integral
    i = int((xx + dmax) / delta)  # x-index
    j = int((yy + dmax) / delta)  # y-index

    #The if-statement places values in the lists
    if ((xx + dmax)/delta) >= 0 and i < num_bins and ((yy + dmax)/delta) >= 0 and j < num_bins:
        count[i,j] += 1
        cur_int[i,j] += current
        Theta[i,j] = theta
        xpt.append(xx)
        ypt.append(yy)

    if rec_prot == plimit: break
    line = data.readline()
################################################################################
 # Distrubtion of the stream of total protons
avg_fluence = rec_prot / (math.pi * radius**2)

# Distrubtion of the stream of protons in the image
im_fluence = count.sum() / (4 * dmax**2)

for i in range(num_bins):
    for j in range(num_bins):
        if count[i,j] != 0:
            cur_int[i,j] /= count[i,j]  # average current in cell
            eta[i,j] = avg_fluence * delta**2 / ( 1.0 - 0.5*cur_int[i,j] )**2

        if count[i,j] >= Cmin:
            Fluct[i,j] = 2.0 * ( 1.0 - math.sqrt(avg_fluence * delta**2/count[i,j]) )
            err = math.sqrt(avg_fluence * delta**2 ) / count[i,j]
            Diff[i,j] = ( Fluct[i,j] - cur_int[i,j] ) / err
            DiffCts[i,j] = ( count[i,j] - eta[i,j] ) / math.sqrt(count[i,j])


Dpos = Diff[ count>=Cmin ]
cur_pos = cur_int[ count>=Cmin ]
Flpos = Fluct[ count>=Cmin ]
DiffCtspos = DiffCts[ count>=Cmin ]
etapos = eta[ count>=Cmin ]

################################################################################

print "Total protons in sample: %d" % rec_prot
print "Total sample Average Fluence: %12.5E ; Image Average Fluence: %12.5E" % (avg_fluence, im_fluence)

print "Mean current integral: %12.5E ; Std. Dev. current integral: %12.5E" % (cur_int.mean(), cur_pos.std())
print "Max current integral: %12.5E ; Min current integral: %12.5E" % (cur_pos.max(), cur_pos.min())

print "Mean Fluct.: %12.5E ; Std. Dev. Fluct.: %12.5E" % (Flpos.mean(), Flpos.std())
print "Max Fluct: %12.5E ; Min Fluct: %12.5E" % (Flpos.max(), Flpos.min())

print "Mean difference: %12.5E ; Std. Dev. Difference: %12.5E" % (Dpos.mean(), Dpos.std())
print "Max difference: %12.5E ; Min difference: %12.5E" % (Dpos.max(), Dpos.min())

print ("Mean  predicted counts per bin: %12.5E ; Std. Dev. predicted counts per bin: %12.5E" % (etapos.mean(), etapos.std()))
print "Max predicted counts per bin: %d ; Min predicted counts per bin: %d" % (etapos.max(), etapos.min())

print "Mean counts per bin: %12.5E ; Std. Dev. Counts per bin: %12.5E" % (count.mean(), count.std())
print "Max counts per bin: %d ; Min counts per bin: %d" % (count.max(), count.min())
print "Number of bins with zero protons: %d" % (count.size - count[ count>0 ].size)
print "Number of bins with %d or fewer protons: %d" % (Cmin, count.size - count[ count>Cmin ].size)

print "Mean (Count-Pred)/err: %12.5E ; Std. Dev.: %12.5E" % (DiffCtspos.mean(), DiffCtspos.std())
print "Max (Count-Pred)/err: %12.5E ; Min (Count-Pred)/err: %12.5E" % (DiffCtspos.max(), DiffCtspos.min())


print "Total protons in image: %d ; Total predicted: %d" % (count.sum(), eta.sum())
#########################Printing Data##########################################

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


for i in range(num_bins):
    for j in range(num_bins):

        xx = ru.idx2vec((i,j))
        line = ("%3d  %3d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  " +\
               "%12.5E  %12.5E\n") %\
               (i, j, xx[0], xx[1], \
                Fluct[i,j], cur_int[i,j], Diff[i,j], count[i,j], eta[i,j], DiffCts[i,j] )

        out.write(line)
out.close()

######################Creating Plots############################################
stretch = 13.0 / 9.0 # scaling factor

fig  = plt.figure()
fig.set_figwidth(18.0 * stretch)
fig.set_figheight(12.0)
######################Radiography.png##########################################
''' The plots 'fluence contrast', 'current projection function',
'Fluence Contrast Difference/Noise','Counts/Bin','Predicted Counts/Bin',
and 'Counts Difference/Noise' '''
'''
# Fluence Contrast
ax = fig.add_subplot(2,3,1)
p = ax.pcolormesh(x,y, Fluct, cmap=cm.afmhot, vmin=min(Fluct.min(), cur_int.min()),
                  vmax=max(Fluct.max(), cur_int.max()))
ax.set_ylabel("Y (cm)")
plt.colorbar(p)
ax.set_title("Fluence Contrast")

# Current Projection Function
ax = fig.add_subplot(2,3,2)
p = ax.pcolormesh(x,y, cur_int, cmap=cm.afmhot, vmin=min(Fluct.min(), cur_int.min()),
                  vmax=max(Fluct.max(), cur_int.max()))
plt.colorbar(p)
ax.set_title("Current Projection Function")

#Fluence Contrast Difference/Noise
vm = max(abs(Diff.min()), abs(Diff.max()))
ax = fig.add_subplot(2,3,3)
p = ax.pcolormesh(x,y, Diff, cmap=cm.RdYlGn, vmax=vm, vmin=-vm)
plt.colorbar(p)
ax.set_title("Fluence Contrast Difference/Noise")
'''
# Counts/Bin
ax = fig.add_subplot(1,1,1)
p = ax.pcolormesh(x,y, count, cmap=cm.afmhot, vmin=count.min(),
                  vmax= count.max())
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
plt.colorbar(p)
ax.set_title("Counts/Bin")
fig.savefig("radiography.png", format='png')
'''
# Predicted Counts/Bin
ax = fig.add_subplot(2,3,5)
p = ax.pcolormesh(x,y, eta, cmap=cm.afmhot, vmin=min(eta.min(),
                  count.min()), vmax=max(eta.max(), count.max()))
ax.set_xlabel("X (cm)")
plt.colorbar(p)
ax.set_title("Predicted Counts/Bin")

# Counts Difference/Noise
vm = max(abs(DiffCts.min()), abs(DiffCts.max()))
ax = fig.add_subplot(2,3,6)
p = ax.pcolormesh(x,y, DiffCts, cmap=cm.RdYlGn, vmax=vm, vmin=-vm)
ax.set_xlabel("X (cm)")
plt.colorbar(p)
ax.set_title("Counts Difference/Noise")


fig.savefig("radiography.png", format='png')

######################Radiography_points.png####################################
fig.clf()
# Proton Locations plot
fig.set_figwidth(6.0)
fig.set_figheight(6.0)
ax = fig.add_subplot(1,1,1)
ax.plot(xpt, ypt, '.', markersize = 0.01 * math.sqrt(1.0E+06/6118107))
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
ax.set_title("Proton Locations")
fig.savefig("radiography_points.png", format='png')


######################Radiography_angdev.png##########################################
fig.clf()
# Angular Deviation
fig.set_figwidth(6.0/0.85)
fig.set_figheight(6.0)
ax = fig.add_subplot(1,1,1)
p = ax.pcolormesh(x,y, Theta, cmap=cm.RdYlGn)
ax.set_xlabel("X (cm)")
ax.set_ylabel("Y (cm)")
ax.set_title("Angular Deviation")
plt.colorbar(p)

fig.savefig("radiography_angdev.png", format='png')
'''
