# -*- coding: utf-8 -*-⏊
import rad_ut as ru
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

x = np.array(range(1,11))
y = np.array(range(20,26))
X, Y = np.meshgrid(x,y)



    print ("Constructing Reconstructed B_perp Projection Plot")
    #BrMag = magnetic_field(Br)
    X,Y = ru.position(flux_ref, bin_um)
    # Intiating Plot
    fig  = plt.figure()
    fig.set_figwidth(7.7)
    fig.set_figheight(6.0)
    # Reconstructed Magnetic Field Plot
