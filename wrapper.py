import sys
sys.path.append(r"/Users/asbogale/Documents/Work/FLASHLAB/pradreader")
from pradreader import reader
import fun_image as image
import fun_Bplot2 as plot
import fun_reconstruct as recon

s2r_cm,s2d_cm,Ep_MeV,flux2D,bin_num = reader.reader(r"/Users/asbogale/Documents/Work/FLASHLAB/ALEMAYEHU/01 - 15 MeV/cshoc_ProtonDetectorFile01_2.201E-09","flash4")

print "flux2D:",d,"\n"
print "s2r_cm:",a,"\n"
print "s2d_cm:",b,"\n"
print "Ep_MeV:",c,"\n"
print "bin_um:",e,"\n"

def wrapper(fn, rtype):
