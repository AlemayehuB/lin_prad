

'''
Acts as a wrapper that runs the recons in this project
'''
import sys

from pradreader import reader
import Bplot2 as plot
import rad_ut as ru
import alogrithm as alog
import path
import image

import  numpy as np
def L2(bin_um, s2r_cm, s2d_cm, BperpR, BperpS):
    '''
    Determines the relative L2 between two arrays

    Parameters
    ----------
    bin_um (float): Length of the side of a bin, in cm
    s2r_cm (float):  Distance from the proton source to the interaction region, in cm
    s2d_cm (float): Distance from the proton source to the detector, in cm
    BperpR (2D array of (x,y)): Reconstructed B Field per (x,y)
    BperpS (2D array of (x,y)): Path integrated B Field per (x,y)

    Returns
    -------
    Constructs a print statement for the user
    '''
    delta_i = (bin_um/10000.0) * s2r_cm / s2d_cm
    EB = BperpS.flatten().dot(BperpS.flatten()) * delta_i**2
    EBR = BperpR.flatten().dot(BperpR.flatten()) * delta_i**2
    L2B = ru.fnorm(BperpR-BperpS)/ru.fnorm(BperpS)

    buf = "EB = %12.5E ; EB_Reconstructed = %12.5E" % (EB, EBR)
    print "...done.  " + buf
    print "Relative L2 norm of (reconstructed B - actual B) = %12.5E" % L2B

def prad_wrap():
    '''
    Wrapper Function for Command line tool that runs the restruction algorithm and
    plotting mechanism

    Parameters
    ----------
    sys.argv[1] (string): The full filename (including path) e.g. "/home/myouts/prad_input"
    sys.argv[2] (string): file type e.g. carlo, flash4, mitcsv

    Returns
    -------
    files (string): png file that contains the reconstructed and/or path integrated
                    magnetic field plot
    '''

    print "STARTING RECONSTRUCTION AND PLOTTING..."
    s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(sys.argv[1],sys.argv[2])
    flux = flux.T
    flux_ref = flux_ref.T

    # Magnetic Field Alogrithm
    print "Calculating Magnetic Perpendicular Field..."
    Bperp = np.zeros((flux.shape[0],flux.shape[0],2))
    if sys.argv[2] == "carlo":
        Bperp, J, avg_fluence, im_fluence = path.mag_parse(r"%s" % sys.argv[1], bin_um)

    BperpR,BperpS = alog.B_recon(flux, flux_ref, Bperp, s2d_cm, s2r_cm, bin_um, Ep_MeV)
    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.B_plot(BperpR, flux_ref, bin_um, sys.argv[2],"Reconstructed")

    # Genereates the Log True B perpendicular Projection,B_true.png
    if sys.argv[2] == "carlo":
        plot.B_plot(BperpS, flux_ref, bin_um, sys.argv[2], "True")

        #Relative L2 between Reconstructed and Path integrated magnetic field
        L2(bin_um, s2r_cm, s2d_cm, BperpR, BperpS)

if __name__=="__main__":
    prad_wrap()
