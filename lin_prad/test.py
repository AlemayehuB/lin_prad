

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

s2r_cm,s2d_cm,Ep_MeV,flux,flux_ref,bin_um = reader.reader(sys.argv[1],sys.argv[2])
flux = flux.T
flux_ref = flux_ref.T
image.hist2D_plot(flux, bin_um, sys.argv[2], "counts")
