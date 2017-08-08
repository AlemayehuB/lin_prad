import math
import sys

import fun_rad_ut as ru
import fun_reconstruct as fr
from constants import M_PROTON_G, ESU, C, V_PER_E, MARG


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Floor counts per bin
Cmin = 10.0

def fluct(flux, num_bins,rap, rs, ri, tot_prot):
    '''
    Calculates the fluence contrast per bin and is put into 2D array

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    num_bins (int): Float, size of the square edge lengths with which to divide the detector for binning
    rap (float): Aperature of the cone that is collimated to screen
    rs (float): Lenght from implosion to screen
    ri (float): Length from implosion to interaction region
    tot_prot (int): number of protons shot at the screen

    Returns
    -------
    Fluct (2D array): fluecne contrast per bin
    '''
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
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    num_bins (int): Float, size of the square edge lengths with which to divide the detector for binning
    rap (float): Aperature of the cone that is collimated to screen
    rs (float): Lenght from implosion to screen
    ri (float): Length from implosion to interaction region
    tot_prot (int): number of protons shot at the screen

    Returns
    -------
    Fluence_Contrast.png
    '''
    print "Constructing Fluence Contrast Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }
    Fluct = fluct(flux,num_bins,rap, rs, ri, tot_prot)
    x = ru.position(num_bins, rap, rs, ri)[0]
    y = ru.position(num_bins, rap, rs, ri)[1]
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, Fluct, cmap=cm.afmhot,
                        vmin=Fluct.min(),vmax=Fluct.max())
    ax.set_ylabel("Y (cm)", fontdict=font)
    ax.set_xlabel("X (cm)", fontdict=font)
    plt.colorbar(p)
    ax.set_title("Fluence Contrast", fontdict=font)
    fig.savefig("Fluence_Contrast.png", format='png')


def flux_plot(flux,num_bins,rap, rs, ri):
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    num_bins (int): Float, size of the square edge lengths with which to divide the detector for binning
    rap (float): Aperature of the cone that is collimated to screen
    rs (float) : Lenght from implosion to screen
    ri (float): Length from implosion to interaction region

    Returns
    -------
    Flux.png
    '''
    print "Constructing Counts/Bin Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }
    x = ru.position(num_bins, rap, rs, ri)[0]
    y = ru.position(num_bins, rap, rs, ri)[1]
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    # Counts/Bin
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, flux, cmap=cm.afmhot, vmin=flux.min(),
                      vmax= flux.max())
    ax.set_xlabel("X (cm)",fontdict=font)
    ax.set_ylabel("Y (cm)",fontdict=font)
    plt.colorbar(p)
    ax.set_title("Counts/Bin",fontdict=font)
    fig.savefig("Flux.png", format='png')
