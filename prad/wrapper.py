#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Acts as a wrapper that links all the separtae modules in this project
'''
import sys
sys.path.append(r"/Users/asbogale/Documents/Work/FLASHLAB/pradreader")
from pradreader import reader
import image as image
import Bplot2 as plot
import reconstruct as fr
import rad_ut as ru

def prad_wrap():
    # First Parameter: Path name of the file
    # Second Parameter: The type of experimental output
    s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(r"%s" % sys.argv[1], "%s" % sys.argv[2])

    # Genereates the Fluence Contrast Plot, fluence_contrast.png
    image.fluct_plot(flux, flux_ref, bin_um, sys.argv[2])

    # Generates the Flux Plot, flux.png
    image.flux_plot(flux, bin_um, sys.argv[2])

    # Reconstructed Magnetic Field Alogrithm
    Br = fr.B_recon(flux, flux_ref, s2d_cm, s2r_cm, bin_um, Ep_MeV)

    #  Genereates the Log Reconstructed B perpendicular Projection,B_recon.png
    plot.BR_plot(Br, flux_ref, bin_um, sys.argv[2])

if __name__=='__main__':
    prad_wrap()
