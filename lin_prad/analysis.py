

'''
Acts as a wrapper that runs the anlysis of the input and plots the data in this project
'''
import sys
import math

from pradreader import reader
import alogrithm as alog
import Bplot2 as plot
import image
import path
import rad_ut as ru

import numpy as np
import pandas as pd


def prad_wrap():
    # Floor counts per bin
    flux_min = 10.0
    x = pd.read_csv(sys.argv[1], header=None,delim_whitespace=True, comment='#',usecols=[0])
    num_prot = len(x)

    print "STARTING ANALYSIS AND PLOTTING..."
    # First Parameter: Path name of the file
    # Second Parameter: The type of experimental output
    s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(sys.argv[1],sys.argv[2])
    flux = flux.T
    flux_ref = flux_ref.T

    # Protons per bin 2D Histogram
    image.hist2D_plot(flux, bin_um, sys.argv[2],"Counts per Bin")
    print "Mean counts per bin: %12.5E ; Std. Dev. Counts per bin: %12.5E" % (flux.mean(), flux.std())
    print "Max counts per bin: %d ; Min counts per bin: %d" % (flux.max(), flux.min())
    print "Number of bins with zero protons: %d" % (flux.size - flux[ flux>0 ].size)
    print "Number of bins with %d or fewer protons: %d" % (flux_min, flux.size - flux[ flux>flux_min ].size)

    # Fluence Distrubtion of protons at the screen 2D Histogram
    Src, fluc = alog.steady_state(flux, flux_ref)
    image.hist2D_plot(fluc, bin_um, sys.argv[2], "Fluence Contrast")
    Flpos = fluc[flux >= flux_min]
    print "Mean Fluct.: %12.5E ; Std. Dev. Fluct.: %12.5E" % (Flpos.mean(), Flpos.std())
    print "Max Fluct: %12.5E ; Min Fluct: %12.5E" % (Flpos.max(), Flpos.min())

    if sys.argv[2] == "carlo":
        x = pd.read_csv(sys.argv[1], header=None,delim_whitespace=True, comment='#',usecols=[0])
        num_prot = len(x)

        Bperp, J, avg_fluence, im_fluence = path.mag_parse(r"%s" % sys.argv[1], bin_um)

        # Current Projection Function 2D Histogram
        image.hist2D_plot(J, bin_um, sys.argv[2],"Current Projection Function")
        Jpos = J[flux >= flux_min]

        # Fluence Cotrast Noise in a 2D Histogram
        err = np.divide(math.sqrt(avg_fluence * (bin_um/10000)**2 ), flux)
        Diff = np.divide(np.subtract(fluc, J), err)
        image.err2D_plot(Diff, bin_um, sys.argv[2],"Fluence Contrast")
        Dpos = Diff[ flux>=flux_min ]

        # Predicted Protons per bin 2D Histogram
        eta = np.divide(avg_fluence * (bin_um/10000.0)**2, np.power((np.subtract(1.0, np.multiply(0.5, J))), 2))
        image.hist2D_plot(eta, bin_um, sys.argv[2], "Predicted Counts per Bin")
        etapos = eta[ flux>=flux_min ]

        # Counts/Bin Differnce/Noise in 2D Histogram
        DiffCts = np.divide(np.subtract(flux, eta), np.sqrt(flux))
        DiffCtspos = DiffCts[flux>=flux_min]
        image.err2D_plot(DiffCts, bin_um, sys.argv[2], "Counts per Bin")

        print "Mean J: %12.5E ; Std. Dev. J: %12.5E" % (Jpos.mean(), Jpos.std())
        print "Max J: %12.5E ; Min J: %12.5E" % (Jpos.max(), Jpos.min())
        print "Mean difference: %12.5E ; Std. Dev. Difference: %12.5E" % (Dpos.mean(), Dpos.std())
        print "Max difference: %12.5E ; Min difference: %12.5E" % (Dpos.max(), Dpos.min())
        print "Mean  predicted counts per bin: %12.5E ; Std. Dev. predicted counts per bin: %12.5E" % (etapos.mean(), etapos.std())
        print "Max predicted counts per bin: %d ; Min predicted counts per bin: %d" % (etapos.max(), etapos.min())
        print "Mean (Count-Pred)/err: %12.5E ; Std. Dev.: %12.5E" % (DiffCtspos.mean(), DiffCtspos.std())
        print "Max (Count-Pred)/err: %12.5E ; Min (Count-Pred)/err: %12.5E" % (DiffCtspos.max(), DiffCtspos.min())
        print "Total sample Average Fluence: %12.5E ; Image Average Fluence: %12.5E" % (avg_fluence, im_fluence)
        print "Total protons in sample: %d" % num_prot

if __name__=="__main__":
    prad_wrap()
