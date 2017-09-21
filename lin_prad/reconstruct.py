

'''
Acts as a wrapper that runs the recons in this project
'''
import sys
import os.path

from pradreader import reader
import Bplot2 as plot
import rad_ut as ru
import alogrithm as alog
import path
import image

import  numpy as np
import argparse as ap

def get_input_data():
    '''
    Command line options and variables
    '''
    avail_input_formats = ['carlo', 'flash4', 'mitcsv']
    parser = ap.ArgumentParser(
                description="This script is used to reconstruct the magnetic field of "
                             "Proton Radiography experiment")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
    parser.add_argument("filename", type=str,
                help="The filename including the path")
    parser.add_argument("filetype", type=str,
                help="The file format",choices=avail_input_formats)
    parser.add_argument("bin_um", type=float,
                help="Length of the side of a bin, in cm")
    parser.add_argument("--tol", default=1.0E-04, type=float,
                help="The Gauss-Seidel tolerance. DEFAULT:1.0E-04")
    parser.add_argument("--iter", default=4000, type=int,
                help="The number of Gauss-Seidel iterations. DEFAULT:4000")

    args = parser.parse_args()

    return args

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
    # Input variables and options
    args = get_input_data()
    fn = args.filename
    rtype = args.filetype
    bin_um = args.bin_um
    tol_iter = args.tol
    max_iter= args.iter
    #############################
    print "STARTING RECONSTRUCTION AND PLOTTING..."
    fname = reader.reader(fn, rtype, 0, 0, 0, 0, bin_um)
    flux,flux_ref,s2r_cm,s2d_cm,Ep_MeV,bin_um = path.parse_im(fname)
    flux = flux.T
    flux_ref = flux_ref.T

    # Magnetic Field Alogrithm
    print "Calculating Magnetic Perpendicular Field..."
    Bperp = np.zeros((flux.shape[0],flux.shape[0],2))

    if rtype == "carlo":
        Bperp, J, avg_fluence, im_fluence = path.mag_parse(fn, bin_um)

    BperpR,BperpS = alog.B_recon(flux, flux_ref, Bperp, s2d_cm, s2r_cm, bin_um, Ep_MeV, tol_iter, max_iter)

    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.B_plot(BperpR, flux_ref, bin_um, rtype,"Reconstructed")

    # Genereates the Log True B perpendicular Projection,B_true.png
    if sys.argv[2] == "carlo":
        plot.B_plot(BperpS, flux_ref, bin_um, rtype, "True")

        #Relative L2 between Reconstructed and Path integrated magnetic field
        L2(bin_um, s2r_cm, s2d_cm, BperpR, BperpS)

if __name__=="__main__":
    prad_wrap()
