import math
import sys

import rad_ut as ru
import fun_reconstruct as fr

import numpy as np

def magnetic_field(flux, rs, ri, rap, tot_prot, num_bins, Tkin):
    '''
    The goal is to calculate the magnetic field per bin
    Parameters
    ----------
    flux (2D array):
        Number of protons per bin
    rs (float) :
        Lenght from implosion to screen
    ri (float):
        Length from implosion to interaction region
    rap (float):
        Aperature of the cone that is collimated to screen
    tot_prot (float):
        Number of protons from the original capsule impolsion
    num_bins (int):
        Bin per dimenison
    Tkin (float):
        Kinetic Energy
    Returns
    -------
    BrMag (2D array): B Field per bin
    '''
    BrMag = np.zeros((num_bins,num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            Br = fr.B_Recon(flux, rs, ri, rap, tot_prot, num_bins, Tkin)
            BrMag = 0.5*math.log10(Br[i,j,0]**2 + Br[i,j,1]**2)
    return BrMag

def position(num_bins):
    '''
    Obtain the positon of each bin
    Parameters
    ----------
    num_bins (int):
        Bin per dimenison
    Returns
    -------
    X (2D array): The x position
    Y (2D array): The y position
    '''
    X = np.zeros((nbins,nbins))
    Y = np.zeros((nbins,nbins))
    for i in range(num_bins):
        for j in range(num_bins):
            x = ru.idx2vec((i,j))
            X[i,j] = x[0]
            Y[i,j] = x[1]
    return (X,Y)

def BR_plot(flux, rs, ri, rap, tot_prot, num_bins, Tkin):
    '''
    Genereates the Log Reconstructed B perpendicularn Projection
    Parameters
    ----------
    flux (2D array):
        Number of protons per bin
    rs (float) :
        Lenght from implosion to screen
    ri (float):
        Length from implosion to interaction region
    rap (float):
        Aperature of the cone that is collimated to screen
    tot_prot (float):
        Number of protons from the original capsule impolsion
    num_bins (int):
        Bin per dimenison
    Tkin (float):
        Kinetic Energy
    Returns
    -------
    B_Reconstructed.png
    '''
    BR = fr.B_Recon(flux, rs, ri, rap, tot_prot, num_bins, Tkin)
    BrMag = magnetic_field(flux, rs, ri, rap, tot_prot, num_bins, Tkin)
    X,Y = position(128)
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

x = fr.flux_image(sys.argv[1],128)
BR_plot(x, 1.00000E+02, 1.00000E+01, 2.00000E-01, 10000000, 128, 14.7)
