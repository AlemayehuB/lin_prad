import math
import sys

import fun_rad_ut as ru
import fun_reconstruct as fr
from constants import M_PROTON_G, ESU, C, V_PER_E


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Floor counts per bin
Cmin = 10.0

def fluct_plot(flux, flux_ref, rs, ri, bin_um):
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    rs (float): Distance from the proton source to the detector, in cm
    ri (float): Distance from the proton source to the interaction region, in cm
    bin_um (float): Length of the side of a bin, in cm
    Tkin (float): Kinetic Energy

    Returns
    -------
    Fluence_Contrast.png (image): Fluence Contrast Plot
    '''
    print "Constructing Fluence Contrast Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }

    Fluct =  fr.steady_state(flux, flux_ref, rs, ri)[1]
    #print "Fluct:",Fluct
    x,y = ru.position(flux, bin_um)
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x, y, Fluct, cmap = cm.afmhot,
                        vmin = Fluct.min(), vmax = Fluct.max())
    ax.set_ylabel("Y (cm)", fontdict=font)
    ax.set_xlabel("X (cm)", fontdict=font)
    plt.colorbar(p)
    ax.set_title("Fluence Contrast", fontdict=font)
    fig.savefig("Fluence_Contrast.png", format='png')


def flux_plot(flux,bin_um):
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    rs (float): Distance from the proton source to the detector, in cm
    ri (float): Distance from the proton source to the interaction region, in cm
    bin_um (float): Length of the side of a bin, in cm

    Returns
    -------
    Flux.png (image): Flux Plot
    '''
    print "Constructing Counts/Bin Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }
    x,y = ru.position(flux, bin_um)
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    # Counts/Bin
    ax = fig.add_subplot(1,1,1)
    #print "Flux:",flux
    p = ax.pcolormesh(x,y, flux, cmap=cm.afmhot, vmin=flux.min(),
                      vmax= flux.max())
    ax.set_xlabel("X (cm)",fontdict=font)
    ax.set_ylabel("Y (cm)",fontdict=font)
    plt.colorbar(p)
    ax.set_title("Counts/Bin",fontdict=font)
    fig.savefig("Flux.png", format='png')

# a = fr.flux_image(sys.argv[1],135)
# flux_plot(a,650)
# b = np.zeros((135,135))
# for i in range(135):
#     for j in range(135):
#         b[i,j] = a.mean()
#
# fluct_plot(a, b, 100, 10, 650)
