'''
Many of the variables are in reference to the paper,
Inferring Morphology and Strength of Magnetic Fields From Proton Radiographs,
This script will reconstruct magnetic fields given the data
from a Proton Radiography experiment.
'''
import sys
import math
from re import match

import fun_rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E, MARG

import numpy as np

# Gauss-Seidel iteration tolerance
TOL_ITER = 1.0E-4

# Maximum number of Gauss-Seidel iterations
MAX_ITER = 4000


def b_field(rs, ri, Tkin):
    '''
    Calculates the Uniform Magnetic field.

    Parameters
    ----------
    rs (float): Lenght from implosion to screen
    ri (float): Length from implosion to interaction region
    Tkin (float): Kinetic energy

    Returns
    -------
    Bconst (float):Calculates the B, uniform magnetic field strength
    '''
    v = math.sqrt(2 * (Tkin*V_PER_E) / M_PROTON_G) # Velocity of Proton

    Bconst = M_PROTON_G * C * v / (ESU * (rs-ri)) #Uniform field Strength

    return Bconst


def steady_state(flux, num_bins, rap, rs, ri, tot_prot):
    '''
    The goal is the obtain the steady-state diffusion Equation

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    num_bins (int): Float, size of the square edge lengths with which to divide the detector for binning
    rap (float): Aperature of the cone that is collimated to screen
    rs (float) : Lenght from implosion to screen
    ri (float): Length from implosion to interaction region
    tot_prot (float): Number of protons from the original capsule impolsion

    Returns
    -------
    Lam (2D array): fluence contrast
    Src(2D array): Source term from multiplying the fluence contrast and exp(fluence contrast)
    '''
    # radius of undeflected image of aperture at screen
    radius = rap * rs / ri
    # Distrubtion of the stream of protons
    avg_fluence = tot_prot/(math.pi * radius**2)

    ru.dmax = MARG * radius / math.sqrt(2.0)
    ru.delta = 2.0 * ru.dmax / num_bins
    Lam = np.ones((num_bins,num_bins))
    # Obtaining the fluence contrast
    a = np.multiply(avg_fluence,np.divide((ru.delta**2),flux)) # variable to break the equation into smaller parts
    Lam = np.multiply(2,np.subtract(1,np.sqrt(a))) # setting the fluence contrast
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


def B_Recon(flux, num_bins, rap, rs, ri, tot_prot, Tkin):
    '''
    Produces a reconstructed magnetic field

    Parameters
    ----------
    flux (2D array): Number of protons per bin
    rap (float): Aperature of the cone that is collimated to screen
    rs (float) : Lenght from implosion to screen
    ri (float): Length from implosion to interaction region
    tot_prot (float): Number of protons from the original capsule impolsion
    num_bins (int): Float, size of the square edge lengths with which to divide the detector for binning
    Tkin (float): Kinetic Energy

    Returns
    -------
    B_R (2D array of (x,y)): Reconstructed Magnetic Field
    '''
    print ("Min Pixel Count: %12.5E\nMax Pixel Count: %12.5E\n"
           "Mean Pixel Count: %12.5E"
           %(flux.min(),flux.max(),flux.mean()))
    # RHS of the Steady-State Diffusion Equation and Fluence Contrast
    Src,Lam = steady_state(flux, num_bins, rap, rs, ri, tot_prot)
    # The real component after Lam is transformed then convolved and then inversely transformed
    phi = ru.solve_poisson(Lam)
    # Iterate to solution
    GS = ru.Gauss_Seidel(phi, np.exp(Lam), D, O, Src, talk=20, tol=TOL_ITER, maxiter=MAX_ITER)
    phi *= ru.delta**2
    # Uniform B Field Strength
    Bconst = b_field(rs, ri, Tkin)
    # Reconstructed perpendicular B Fields
    B_R = np.zeros((num_bins, num_bins,2))
    # Lateral motion of proton
    deltaX = np.zeros((num_bins, num_bins,2))

    for i in range(num_bins):
        for j in range(num_bins):
            #Reconstructed Data
            deltaX[i,j] = -ru.gradient(phi, (i,j))
            B_R[i,j,0] = Bconst * deltaX[i,j,1]
            B_R[i,j,1] = -Bconst * deltaX[i,j,0]

    return B_R


def flux_image(filename, num_bins):
    #Data file descriptor
    data = open(filename, 'r')
    plimit =- 1
    line = data.readline()
    while not match('# Columns:', line): #Sets the variables

            if match('# Tkin:', line):
                Tkin = float(line.split()[2]) #Kinetic energy

            elif match('# rs:', line):
                rs = float(line.split()[2]) # Length from implosion to screen

            elif match('# ri:', line):
                ri = float(line.split()[2]) # Length from implosion to interaction region

            elif match('# raperture:', line):
                rap = float(line.split()[2]) # aperture

            line = data.readline()
    # The loop reads each line of the input until the data is reached.
    while match('#', line): line = data.readline()
    # Pixel Counts
    count = np.zeros((num_bins, num_bins))
    rec_prot =  0
    radius = rap * rs / ri  # radius of undeflected image of aperture at screen
    ru.dmax = MARG * radius / math.sqrt(2.0)  # variable in rad_ut
    ru.delta = 2.0 * ru.dmax / num_bins  # variable in rad_ut
    while line:
        rec_prot += 1
        line = line.split()
        x_loc = float(line[3]) # Final X location at screen (cm)
        y_loc = float(line[4]) # Final Y location at screen (cm)
        i,j = ru.vec2idx((x_loc,y_loc))

        #The if-statement places values in the lists
        if (x_loc + ru.dmax)/ru.delta >= 0 and i < num_bins and (y_loc + ru.dmax)/ru.delta >= 0 and j < num_bins:
            count[i,j] += 1

        if rec_prot == plimit: break
        line = data.readline()
    return count
