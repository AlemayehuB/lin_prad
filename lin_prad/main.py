

'''
Acts as a wrapper that links all the separtae modules in this project
'''
import sys
from pradreader import reader
import image as image
import Bplot2 as plot
import reconstruct as fr
import rad_ut as ru

def prad_wrap():
    print "STARTING RECONSTRUCTION and ANALYSIS.."
    # First Parameter: Path name of the file
    # Second Parameter: The type of experimental output
    s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(r"%s" % sys.argv[1], "%s" % sys.argv[2])

    # Genereates the Fluence Contrast Plot, fluence_contrast.png
    image.fluct_plot(flux, flux_ref, bin_um, sys.argv[2])

    # Generates the Flux Plot, flux.png
    image.flux_plot(flux, bin_um, sys.argv[2])

    # Reconstructed Magnetic Field Alogrithm
    print "Calculating Reconstructed Magnetic Perpendicular Field"
    tol_iter = float(raw_input(r"Desired Gauss-Seidel tolerance [default 1.0E-04]:  ") or 1.0E-04)
    max_iter = int(raw_input(r"Desired number of Gauss-Seidel iterations [default 4000]:  ") or 4000)
    Br = fr.B_recon(flux, flux_ref, s2d_cm, s2r_cm, bin_um, Ep_MeV, tol_iter, max_iter)

    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.BR_plot(Br, flux_ref, bin_um, sys.argv[2])

if __name__=="__main__":
    prad_wrap()
