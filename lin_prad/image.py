
'''
Provides tools to analyze the inputted data from a proton radiogtaphy experiment
'''

import math
import sys

import rad_ut as ru
import reconstruct as fr
from constants import M_PROTON_G, ESU, C, V_PER_E


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Floor counts per bin
Cmin = 10.0

def fluct_plot(flux, flux_ref, bin_um, type):
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    bin_um (float): Length of the side of a bin, in cm

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

    fluct =  fr.steady_state(flux, flux_ref)[1]
    x,y = ru.position(flux, bin_um)
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)
    ax = fig.add_subplot(1,1,1)
    print x.shape
    print y.shape
    print fluct.shape
    p = ax.pcolormesh(x, y, fluct, cmap = cm.afmhot,
                        vmin = fluct.min(), vmax = fluct.max())
    ax.set_ylabel("Y (cm)", fontdict=font)
    ax.set_xlabel("X (cm)", fontdict=font)
    plt.colorbar(p)
    if type == 'carlo':
        x = "Carlo"
    elif type == 'flash4':
        x = "Flash"
    elif type == 'mitcsv':
        x = 'MITCSV'
    ax.set_title(x + ": Fluence Contrast", fontdict=font)
    fig.savefig("fluence_contrast.png", format='png')


def flux_plot(flux, bin_um, type):
    '''
    Genereates the a plot that shows the number of
    protons per bin

    Parameters
    ----------
    flux (2D array): Number of protons per bin
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
    p = ax.pcolormesh(x,y, flux, cmap=cm.afmhot, vmin=flux.min(),
                      vmax= flux.max())
    ax.set_xlabel("X (cm)",fontdict=font)
    ax.set_ylabel("Y (cm)",fontdict=font)
    plt.colorbar(p)
    if type == 'carlo':
        x = "Carlo"
    elif type == 'flash4':
        x = "Flash"
    elif type == 'mitcsv':
        x = 'MITCSV'
    ax.set_title(x + ": Counts/Bin",fontdict=font)
    fig.savefig("flux.png", format='png')
