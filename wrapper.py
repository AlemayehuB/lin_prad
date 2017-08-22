#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
sys.path.append(r"/Users/asbogale/Documents/Work/FLASHLAB/pradreader")
from pradreader import reader
import fun_image as image
import fun_Bplot2 as plot
import fun_reconstruct as fr
import fun_rad_ut as ru


if __name__=='__main__':
    # First Parameter: Path name of the file
    # Second Parameter: The type of experimental output
    s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(r"%s" % sys.argv[1], "%s" % sys.argv[2])

    #Genereates the Fluence Contrast Plot, fluence_contrast.png
    image.fluct_plot(flux, flux_ref, s2d_cm, s2r_cm, bin_um)

    # Generates the Flux Plot, flux.png
    image.flux_plot(flux, bin_um)

    # Reconstructed Magnetic Field Alogrithm
    Br = fr.B_recon(flux, flux_ref, s2d_cm, s2r_cm, bin_um, Ep_MeV)

    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.BR_plot(Br, flux_ref, bin_um)
