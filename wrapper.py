import sys
sys.path.append(r"/Users/asbogale/Documents/Work/FLASHLAB/pradreader")
from pradreader import reader
import fun_image as image
import fun_Bplot2 as plot
import fun_reconstruct as fr

ri,rs,Ep_MeV,flux2D,flux2D_ref,bin_um = reader.reader(r"/Users/asbogale/Documents/Work/FLASHLAB/ALEMAYEHU/01 - 15 MeV/cshoc_ProtonDetectorFile01_2.201E-09","flash4")
def wrapper(fn, rtype):
    ri, rs, Ep_MeV, flux2D, flux2D_ref, bin_um = reader.reader(fn, rtype)
    # Reconstruction Field Algorithm
    recon.B_Recon(flux, flux_ref, rs, ri, bin_um, Tkin)

    # Fluence Contrast Plot
    image.fluct_plot(flux, flux_ref, rs, ri, bin_um, Tkin)

    # Flux Plot
    image.flux_plot(flux, flux_ref, rs, ri, bin_um)

    # Genereates the Log Reconstructed B perpendicular Projection
    plot.BR_plot(flux, flux_ref, rs, ri, bin_um, Tkin)


#image.flux_plot(flux2D, bin_um)
#image.fluct_plot(flux2D, flux2D_ref, rs, ri, bin_um)
# fr.B_Recon(flux2D, flux2D_ref, rs, ri, bin_um, Ep_MeV)

plot.BR_plot(flux2D, flux2D_ref, rs, ri, bin_um, Ep_MeV)

#if __name__=='__main__':
