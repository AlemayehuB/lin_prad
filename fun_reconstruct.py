'''
Many of the variables are in reference to the paper,
Inferring Morphology and Strength of Magnetic Fields From Proton Radiographs,
This script will reconstruct magnetic fields given the data
from a Proton Radiography experiment.
'''
from sys import argv
from re import match
import math

import rad_ut as ru
from .constants import M_PROTON_G, ESU, C, V_PER_E
import numpy as np

# margin
MARG = 0.98

# Gauss-Seidel iteration tolerance
TOL_ITER = 1.0E-4

# Maximum number of Gauss-Seidel iterations
MAX_ITER = 4000

def pixel_count()


def b_field(rs, ri, Tkin):
''' Calculates the Uniform Magnetic field.
    Parameters
    ----------
    rs : floar
        Lenght from implosion to screen
    ri : float
        Length from implosion to interaction region
    Tkin : float
        Proton Kinetic Energies
    Returns
    -------
    int
        Calculates the B, uniform field strength
'''
    v = math.sqrt(2 * (Tkin*VperE) / m) # Velocity of Proton

    Bconst = m * c * v / (e * (rs-ri)) #Uniform field Strength

    return Bconst

def steady_state(flux, rs, ri, rap, tot_prot, num_bins):
''' The goal is the obtain the steady-state diffusion Equation
    Parameters
    ----------
    flux : 2D array
        Number of protons per bin
    rs : float
        Lenght from implosion to screen
    ri : float
        Length from implosion to interaction region
    rap : float
        Aperature of the cone that is collimated to screen
    tot_prot: float
        Number of protons from the original capsule impolsion
    num_bins: 2D array
        Bin per dimenison
    Returns
    -------
    returns the source term obtained from multiplying the fluence contrast,Λ, and exp(Λ)
    along with the actual fluence contrast
'''
    # radius of undeflected image of aperture at screen
    radius = rap * rs / ri
    # Distrubtion of the stream of protons
    avg_fluence = tot_prot/(math.pi * radius**2)

    ru.dmax = marg * radius / math.sqrt(2.0)  # variable in rad_ut
    ru.delta = 2.0 * ru.dmax / num_bins  # variable in rad_ut
    # Obtaining the fluence contrast
    Lam = np.zeros((num_bins, num_bins))
    a = np.multiply(avg_fluence,np.divide((ru.delta**2),flux)) # variable to break the equation into smaller parts
    Lam = np.multiply(2,np.subtract(1,np.sqrt(a))) # setting the fluence contrast
    # Obtaining the exponential fluence contrast
    ExpLam = np.exp(Lam)
    # Source Term
    Src = np.multiply(Lam,ExpLam) # RHS of the Steady-State Diffusion Equation
    return (Src,Lam)

def D(i, j, x):
''' Supplemental function used during Gauss-Seidel Iteration
    Parameters
'''
    d = -2.0 * x[i,j] - 0.5 * ( ru.bc_enforce_N(x, i+1,j)
          + ru.bc_enforce_N(x, i-1,j)
          + ru.bc_enforce_N(x, i,j+1)
          + ru.bc_enforce_N(x, i,j-1) )
    return d

def O(i,j, x):
''' Supplemental function used during Gauss-Seidel Iteration
'''
    a = 0.5 * ( ru.bc_enforce_D(x, i+1,j) * (ru.bc_enforce_N(ExpLam, i+1,j) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i-1,j) * (ru.bc_enforce_N(ExpLam, i-1,j) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i,j+1) * (ru.bc_enforce_N(ExpLam, i,j+1) + ExpLam[i,j])
        + ru.bc_enforce_D(x, i,j-1) * (ru.bc_enforce_N(ExpLam, i,j-1) + ExpLam[i,j]))
    return a

def flux_alog(flux,):
'''
'''
