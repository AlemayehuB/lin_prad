

'''
Many of the variables  and equations are in reference to the paper,
Inferring Morphology and Strength of Magnetic Fields From Proton Radiographs,
This script will reconstruct magnetic fields given the data
from a Proton Radiography experiment.
'''
import sys
import math
from re import match

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

import numpy as np

# Gauss-Seidel iteration tolerance
TOL_ITER = 1.0E-4

# Maximum number of Gauss-Seidel iterations
MAX_ITER = 4000


def b_field(s2r_cm, s2d_cm, Ep_MeV):
    '''
    Calculates the Uniform Magnetic field.

    Parameters
    ----------
    s2r_cm (float): Distance from the proton source to the detector, in cm
    s2d_cm (float): Distance from the proton source to the interaction region, in cm
    Ep_MeV (float): Kinetic energy

    Returns
    -------
    Bconst (float): Calculates the B, uniform magnetic field strength
    '''
    v = math.sqrt(2 * (Ep_MeV*V_PER_E) / M_PROTON_G) # Velocity of Proton

    Bconst = M_PROTON_G * C * v / (ESU * (s2r_cm - s2d_cm)) # Uniform B Field Strength

    return Bconst


def steady_state(flux, flux_ref):
    '''
    The goal is the obtain the steady-state diffusion Equation

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    s2r_cm (float): Distance from the proton source to the detector, in cm
    s2d_cm (float): Distance from the proton source to the interaction region, in cm

    Returns
    -------
    Lam (2D array): fluence contrast
    Src (2D array): Source term from multiplying the fluence contrast and exp(fluence contrast)
    '''
    num_bins = flux_ref.shape[0] # num_bins x num_bins
    Lam = np.ones((num_bins,num_bins))
    # Obtaining the fluence contrast from Equation 6
    Lam = np.divide(np.subtract(flux,flux_ref), flux_ref)
    # Obtaining the exponential fluence contrast
    ExpLam = np.exp(Lam)
    # Source Term
    Src = np.multiply(Lam,ExpLam) # RHS of the Steady-State Diffusion Equation

    return (Src,Lam)


def D(i, j, y):
    '''
    Supplemental function used during Gauss-Seidel Iteration
    '''
    d = -2.0 * y[i,j] - 0.5 * ( ru.bc_enforce_N(y, i+1,j)
          + ru.bc_enforce_N(y, i-1,j)
          + ru.bc_enforce_N(y, i,j+1)
          + ru.bc_enforce_N(y, i,j-1) )

    return d


def O(i, j, x, y):
    '''
    Supplemental function used during Gauss-Seidel Iteration
    '''
    a = 0.5 * ( ru.bc_enforce_D(x, i+1,j) * (ru.bc_enforce_N(y, i+1,j) + y[i,j])
        + ru.bc_enforce_D(x, i-1,j) * (ru.bc_enforce_N(y, i-1,j) + y[i,j])
        + ru.bc_enforce_D(x, i,j+1) * (ru.bc_enforce_N(y, i,j+1) + y[i,j])
        + ru.bc_enforce_D(x, i,j-1) * (ru.bc_enforce_N(y, i,j-1) + y[i,j]))

    return a


def B_recon(flux, flux_ref, s2r_cm, s2d_cm, bin_um, Ep_MeV, tol_iter, max_iter):
    '''
    Produces a reconstructed magnetic field

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    s2r_cm (float): Distance from the proton source to the detector, in cm
    s2d_cm (float): Distance from the proton source to the interaction region, in cm
    bin_um (float): Length of the side of a bin, in cm
    Ep_MeV (float): Kinetic Energy

    Returns
    -------
    Br (2D array of (x,y)): Reconstructed Magnetic Field
    '''

    ru.delta = bin_um
    num_bins = flux_ref.shape[0] # num_bins x num_bins
    # RHS of the Steady-State Diffusion Equation and Fluence Contrast
    Src,Lam = steady_state(flux, flux_ref)
    # The real component after Lam is transformed then convolved and then inversely transformed
    phi = ru.solve_poisson(Lam)
    # Iterate to solution
    print "Gauss-Seidel Iteration..."
    GS = ru.Gauss_Seidel(phi, np.exp(Lam), D, O, Src, talk=20, tol=tol_iter, maxiter= max_iter)
    # Multiplying by the area of the bin
    phi *= (ru.delta**2)
    # Uniform B Field Strength
    Bconst = b_field(s2r_cm, s2d_cm, Ep_MeV)
    # Reconstructed perpendicular B Fields
    Br = np.zeros((num_bins, num_bins,2))
    # Lateral motion of proton
    deltaX = np.zeros((num_bins, num_bins,2))

    for i in range(num_bins):
        for j in range(num_bins):
            #Reconstructed Datam
            deltaX[i,j] = -ru.gradient(phi, (i,j))
            Br[i,j,0] = Bconst * deltaX[i,j,1]
            Br[i,j,1] = -Bconst * deltaX[i,j,0]

    return Br
