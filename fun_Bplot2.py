import math
import sys

import fun_rad_ut as ru
import fun_reconstruct as fr
from constants import M_PROTON_G, ESU, C, V_PER_E

import matplotlib.pyplot as plt
import numpy as np

def magnetic_field(flux, flux_ref, rs, ri, bin_um, Tkin):
    '''
    The goal is to calculate the magnetic field per bin

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
    BrMag (2D array): B Field per bin
    '''
    num_bins = flux_ref.shape[0] # num_bins x num_bins
    BrMag = np.zeros((num_bins,num_bins))
    # Reconstructed Magnetic Field
    Br = fr.B_Recon(flux, flux_ref, rs, ri, bin_um, Tkin)
    for i in range(num_bins):
        for j in range(num_bins):
            # Field Strength on a logarithmic scale
            BrMag = 0.5*math.log10(Br[i,j,0]**2 + Br[i,j,1]**2)

    return BrMag


def BR_plot(flux, flux_ref, rs, ri, bin_um, Tkin):
    '''
    Genereates the Log Reconstructed B perpendicular Projection

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
    B_Reconstructed.png (image): Log Reconstructed B perpendicular Projection
    '''
    print (r"Constructing Reconstructed $B_\perp$ Projection Plot")
    BR = fr.B_Recon(flux, flux_ref, rs, ri, bin_um, Tkin)
    BrMag = magnetic_field(flux, flux_ref, rs, ri, bin_um, Tkin)
    X,Y = ru.position(flux_ref, bin_um, rs, ri)
    # Intiating Plot
    fig  = plt.figure()
    fig.set_figwidth(7.7)
    fig.set_figheight(6.0)
    # Reconstructed Magnetic Field Plot
    ax = fig.add_subplot(1,1,1)
    strm = ax.streamplot(X[:,0], Y[0,:], BR[:,:,0].T, BR[:,:,1].T, color=BrMag.T, \
                          linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0)
    fig.colorbar(strm.lines)
    ax.set_title(r"Log Reconstructed $B_\perp$ Projection (G cm)", fontsize=18)
    ax.set_xlabel(r"X (cm)", fontsize=18)
    ax.set_ylabel(r"Y (cm)", fontsize=18)
    ax.tick_params(labelsize='large')

    fig.savefig("B_Reconstructed.png", format='png')
