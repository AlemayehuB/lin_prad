

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

def prad_wrap():
    '''Wrapper Function for command line reconstruction tool'''

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
    for i in range(BperpR.shape[0]):
        for j in range(BperpR.shape[0]):
            print BperpR[i,j],BperpS[i,j]
    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.B_plot(BperpR, flux_ref, bin_um, sys.argv[2],"Reconstructed")

    # Genereates the Log True B perpendicular Projection,B_true.png
    if sys.argv[2] == "carlo":
        plot.B_plot(BperpS, flux_ref, bin_um, sys.argv[2], "True")

if __name__=="__main__":
    prad_wrap()
