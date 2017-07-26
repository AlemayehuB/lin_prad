'''
Many of the variables are in reference to the paper,
Inferring Morphology and Strength of Magnetic Fields From Proton Radiographs,
This script will reconstruct magnetic fields given the data
from a Proton Radiography experiment.
'''
from sys import argv
from re import match
import math

import rad_ut as ru

import numpy as np


#Data file descriptor
data = open(argv[1], 'r')

#Output file descriptor
out = open(argv[2], 'w')

#Proton limiter
plimit =- 1
if len(argv) > 3: plimit = int(argv[3])

#Resolution of plot
num_bins = 128
if len(argv) > 4: num_bins = int(argv[4])

#############Parsing Data File###################################################
out.write("# Input filename: %s\n" % argv[1])
line = data.readline()
while not match('# Columns:', line): #Sets the variables

        if match('# Tkin:', line):
            Tkin = float(line.split()[2]) #Kinetic energy

        elif match('# rs:', line):
            rs = float(line.split()[2]) # Length from implosion to screen

        elif match('# ri:', line):
            ri = float(line.split()[2]) # Length from implosion to interaction region

        elif match('# raperture:', line):
            rap = float(line.split()[2]) # aperture

        out.write(line)
        line = data.readline()

# The loop reads each line of the input until the data is reached.
while match('#', line): line = data.readline()

#############Setting Varibales###########################################
#margin
marg = 0.98

# Gauss-Seidel iteration tolerance
tol_iter = 1.0E-4

# Maximum number of Gauss-Seidel iterations
max_iter = 4000

radius = rap * rs / ri  # radius of undeflected image of aperture at screen
ru.dmax = marg * radius / math.sqrt(2.0)  # variable in rad_ut
ru.delta = 2.0 * ru.dmax / num_bins  # variable in rad_ut
delta_i = ru.delta * ri / rs

e = 4.8032E-10   # Statcoul
m = 1.6726E-24  # mass of proton (g)
c = 2.9979E+10   # Speed of light (cm/s)
VperE= 1.6022E-06 # 1 MeV / 1 erg
v = math.sqrt(2 * (Tkin*VperE) / m) # Velocity of Proton
Bconst = m * c * v / (e * (rs-ri)) #Uniform field Strength

rec_prot=0 #number of recovered protons that hit the collimated square image

# Pixel Counts
count = np.zeros((num_bins, num_bins))
# Current integrals
cur_int = np.zeros((num_bins, num_bins))
# Position-corrected Current integrals
cur_intS = np.zeros((num_bins, num_bins))
# Fluence Contrast
Lam = np.zeros((num_bins, num_bins))
# Position-corrected Fluence Contrast
LamS = np.zeros((num_bins, num_bins))
# B Field Integral
B = np.zeros((num_bins, num_bins, 2))
# Reconstructed perpendicular B Fields
B_R = np.zeros((num_bins, num_bins,2))
#True perpendicular B Fields
B_S = np.zeros((num_bins, num_bins,2))
# Lateral motion of proton
deltaX = np.zeros((num_bins, num_bins,2))
# Reconstructed Lateral motion of proton
deltaXR = np.zeros((num_bins, num_bins,2))
# True Lateral motion of proton
deltaXS = np.zeros((num_bins, num_bins,2))

################################################################################
while line:
    rec_prot += 1
    line = line.split()
    x_i = float(line[0]) # Initial normal vector X-component
    y_i = float(line[1]) # Initial normal vector Y-component
    x_loc = float(line[3]) # Final X location at screen (cm)
    y_loc = float(line[4]) # Final Y location at screen (cm)
    current = float(line[8]) # Current integral
    Bx = float(line[9]) # Bx Integral
    By = float(line[10]) # By Integral
    i,j = ru.vec2idx((x_loc,y_loc))

    #The if-statement places values in the lists
    if (x_loc + ru.dmax)/ru.delta >= 0 and i < num_bins and (y_loc + ru.dmax)/ru.delta >= 0 and j < num_bins:
        count[i,j] += 1
        B[i,j,0] += Bx
        B[i,j,1] += By
        deltaX[i,j,0] += (x_loc - x_i*rs)
        deltaX[i,j,1] += (y_loc - y_i*rs)
        cur_int[i,j] += current

    if rec_prot == plimit: break
    line = data.readline()
data.close()

 #Distrubtion of the stream of protons
avg_fluence = rec_prot/(math.pi * radius**2)

################################################################################
buf = "# num_bins = %d ; delta =  %12.5E\n" % (num_bins, ru.delta)
out.write(buf)
print "Reading %s, binning %d x %d..." % (argv[1], num_bins, num_bins)

buf = "Read %d protons, of which %d in square window." % (rec_prot, count.sum())
print buf
out.write("# %s\n" % buf)

buf = ("Min Pixel Count: %12.5E\nMax Pixel Count: %12.5E\n"
       "Mean Pixel Count: %12.5E\nDelta: %12.5E "
       %(count.min(),count.max(),count.mean(),ru.delta))
print buf

#############Matrix Operations#######################################################
try:
    ##########B Field################
    B.shape=(num_bins*num_bins,2)
    count.shape=(num_bins*num_bins,1)
    B = np.divide(B,count)
    B.shape=(num_bins,num_bins,2)

    ##########deltaX##################
    deltaX.shape=(num_bins*num_bins,2)
    deltaX = np.divide(deltaX,count)
    deltaX.shape=(num_bins,num_bins,2)

    #########Current Integral#################
    count.shape=(num_bins, num_bins)
    cur_int = np.divide(cur_int,count)

    #########Lams#####################
    pre_Lam=np.multiply(avg_fluence,np.divide((ru.delta**2),count))
    Lam=np.multiply(2,np.subtract(1,np.sqrt(pre_Lam)))
    ExpLam = np.exp(Lam)
    Src = np.multiply(Lam,ExpLam) #Steady-State Diffusion Equation

except ZeroDivisionError:
    print "Zero pixel, will screw everything up."
    raise ValueError
#############Gauss-Seidel iterations#########################################################
print "Gauss-Seidel Iteration..."

def D(i,j):
    d = -2.0 * ExpLam[i,j] - 0.5 * ( ru.bc_enforce_N(ExpLam, i+1,j)
          + ru.bc_enforce_N(ExpLam, i-1,j)
          + ru.bc_enforce_N(ExpLam, i,j+1)
          + ru.bc_enforce_N(ExpLam, i,j-1) )
    return d

def O(i,j, x):
    a = 0.5 * ( ru.bc_enforce_D(x, i+1,j) * (ru.bc_enforce_N(ExpLam, i+1,j) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i-1,j) * (ru.bc_enforce_N(ExpLam, i-1,j) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i,j+1) * (ru.bc_enforce_N(ExpLam, i,j+1) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i,j-1) * (ru.bc_enforce_N(ExpLam, i,j-1) + ExpLam[i,j]))
    return a

# Initial guess by FFT Poisson solve
phi = ru.solve_poisson(Lam)

# Iterate to solution
GS = ru.Gauss_Seidel(phi, D, O, Src, talk=20, tol=tol_iter, maxiter=max_iter)

phi *= ru.delta**2


for i in range(num_bins):
    for j in range(num_bins):
        #Reconstructed Data
        deltaXR[i,j] = -ru.gradient(phi, (i,j))
        B_R[i,j,0] = Bconst * deltaXR[i,j,1]
        B_R[i,j,1] = -Bconst * deltaXR[i,j,0]

        #True Data
        x = ru.idx2vec((i,j))
        x = x + deltaXR[i,j]
        idx = ru.vec2idx(x)
        ii = idx[0] % num_bins
        jj = idx[1] % num_bins
        B_S[i,j,:] = B[ii, jj, :]
        deltaXS[i,j,:] = deltaX[ii, jj, :]
        cur_intS[i,j] = cur_int[ii, jj]
        LamS[i,j] = Lam[ii, jj]


EB = B_S.flatten().dot(B_S.flatten()) * delta_i**2 # True B Field
EBR = B_R.flatten().dot(B_R.flatten()) * delta_i**2 # Reconstructed B Field
L2B = ru.fnorm(B_R-B_S)/ru.fnorm(B_S) # L2 norm

#############Output Data################################################
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
out.write('# Column 3: x position\n')
out.write('# Column 4: y position\n')
out.write('# Column 5: Bx (True)\n')
out.write('# Column 6: By (True)\n')
out.write('# Column 7: Bx (Reconstructed)\n')
out.write('# Column 8: By (Reconstructed)\n')
out.write('# Column 9: Delta_x (True)\n')
out.write('# Column 10: Delta_y (True)\n')
out.write('# Column 11: Delta_x (Reconstructed)\n')
out.write('# Column 12: Delta_y (Reconstructed)\n')
out.write('# Column 13: Position-corrected Lambda\n')
out.write('# Column 14: Position-corrected current integral\n')

for i in range(num_bins):
    for j in range(num_bins):

        x = ru.idx2vec((i,j))
        line = ("%3d  %3d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  " +\
               "%12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n") %\
               (i, j, x[0], x[1], \
                B_S[i,j,0], B_S[i,j,1], \
                B_R[i,j,0], B_R[i,j,1], \
                deltaXS[i,j,0], deltaXS[i,j,1], \
                deltaXR[i,j,0], deltaXR[i,j,1], \
                LamS[i,j], cur_intS[i,j] )

        out.write(line)
out.close()
