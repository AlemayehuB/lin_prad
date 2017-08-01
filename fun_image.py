import math
import sys

import fun_reconstruct as fr

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
# margin
MARG = 0.98

# Floor counts per bin
Cmin = 10.0

def list(num_bins, rap, rs, ri):
    radius = rap * rs / ri  # radius of the undeflected image on the screen
    dmax = MARG * radius / math.sqrt(2.0)
    delta = 2.0 * dmax / num_bins
    x = np.zeros((num_bins+1, num_bins+1))
    y = np.zeros((num_bins+1, num_bins+1))
    for i in range(num_bins+1):
        xx = -dmax + i*delta
        for j in range(num_bins+1):
            yy = -dmax + j*delta
            x[i,j] = xx
            y[i,j] = yy
    return (x,y)

def fluct(flux, num_bins,rap, rs, ri, tot_prot):
    # Fluence Contrast
    Fluct = np.zeros((num_bins,num_bins))
    radius = rap * rs / ri  # radius of the undeflected image on the screen
    dmax = MARG * radius / math.sqrt(2.0)
    delta = 2.0 * dmax / num_bins
    avg_fluence = tot_prot / (math.pi * radius**2)
    for i in range(num_bins):
        for j in range(num_bins):
            if flux[i,j] >= Cmin:
                Fluct[i,j] = 2.0 * ( 1.0 - math.sqrt(avg_fluence * delta**2/flux[i,j]))
    return Fluct

def fluct_plot(flux, num_bins,rap, rs, ri, tot_prot):
    Fluct = fluct(flux,num_bins,rap, rs, ri, tot_prot)
    x = list(num_bins, rap, rs, ri)[0]
    y = list(num_bins, rap, rs, ri)[1]
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, Fluct, cmap=cm.afmhot,
                        vmin=Fluct.min(),vmax=Fluct.max())
    ax.set_ylabel("Y (cm)")
    ax.set_xlabel("X (cm)")
    plt.colorbar(p)
    ax.set_title("Fluence Contrast")
    print "Making Fluence Contrast Plot"
    fig.savefig("Fluence_Contrast.png", format='png')

def flux_plot(flux,num_bins,rap, rs, ri):
    x = list(num_bins, rap, rs, ri)[0]
    y = list(num_bins, rap, rs, ri)[1]
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    # Counts/Bin
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, flux, cmap=cm.afmhot, vmin=flux.min(),
                      vmax= flux.max())
    ax.set_xlabel("X (cm)")
    ax.set_ylabel("Y (cm)")
    plt.colorbar(p)
    ax.set_title("Counts/Bin")
    fig.savefig("Counts/Bin.png", format='png')
    #fig.savefig("Counts/Bin.png", format='png')

x = fr.flux_image(sys.argv[1],128)
fluct_plot(x,128,  2.00000E-01,1.00000E+02,1.00000E+01, 10000000)
