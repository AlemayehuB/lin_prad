#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Provides function that will plot the reconstructed magnetic field caluctlated from
the alogrithm in recconstruct.py
'''

import math
import sys

import fun_rad_ut as ru
import fun_reconstruct as fr
from constants import M_PROTON_G, ESU, C, V_PER_E

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

def magnetic_field(Br):
    '''
    The goal is to calculate the magnetic field per bin

    Parameters
    ----------
    Br (2D array of (x,y)): B Field per (x,y)

    Returns
    -------
    BrMag (2D array): Log reconstructed B Field per bin
    '''
    num_bins = Br.shape[0] # num_bins x num_bins
    BrMag = np.zeros((num_bins,num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            # Field Strength on a logarithmic scale
            BrMag[i,j]= 0.5*math.log10(Br[i,j,0]**2 + Br[i,j,1]**2)

    return BrMag


def BR_plot(Br, flux_ref, bin_um):
    '''
    Genereates the Log Reconstructed B perpendicular Projection

    Parameters
    ----------
    Br (2D array of (x,y)): B Field per (x,y)
    flux_ref (2D array): Number protons per bin without an interaction region
    bin_um (float): Length of the side of a bin, in cm

    Returns
    -------
    B_recon.png (image): Log Reconstructed B perpendicular Projection
    '''
    print ("Constructing Reconstructed B_perp Projection Plot")
    BrMag = magnetic_field(Br)
    X,Y = ru.position(flux_ref, bin_um)
    # Intiating Plot
    fig  = plt.figure()
    fig.set_figwidth(7.7)
    fig.set_figheight(6.0)
    # Reconstructed Magnetic Field Plot
    ax = fig.add_subplot(1,1,1)
    strm = ax.streamplot(X[:,0], Y[0,:], Br[:,:,0].T, Br[:,:,1].T, color=BrMag.T, \
                          linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0)
    fig.colorbar(strm.lines)
    ax.set_title(r"Log Reconstructed $B_\perp$ Projection (G cm)", fontsize=18)
    ax.set_xlabel(r"X (cm)", fontsize=18)
    ax.set_ylabel(r"Y (cm)", fontsize=18)
    ax.tick_params(labelsize='large')

    fig.savefig("B_recon.png", format='png')
