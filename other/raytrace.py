# Support classes for proton radiographic ray tracing
#
# Carlo Graziani, University of Chicago
#
"""
This module provides classes to support proton radiographic ray-tracing.
"""

import numpy as np
import math
from scipy.stats import uniform
from scipy.integrate import odeint
from scipy.integrate import ode
from mpi4py import MPI
import re
import sys

################################################################################
class pr_magfield(object):

    """
    Subclass this to provide a function returning the RHS of the ODE for
    the ray.

    Constructor invocation:
    pr_magfield(b0, li, Tkin, euler_angles=None, Jint=False, Bperpint=False, switch=False, **kwargs)

    b0:    A typical field value, used to scale the returned field, and to
           set the Larmor radius (Gauss)

    li:    The interaction region length, along the proton ray direction (cm)

    Tkin:  Proton kinetic energy (MeV)

    euler_angles: tuple of 3 rotation angles to characterize the rigid rotation of the
                  field configuration with respect to the fixed proton radiographic
                  frame

    Jint: If true, RHS function also returns integral of current along ray.

    Bperpint: If true, RHS function also returns integral of Bperp along ray.

    dBperpint: If true, RHS function also returns integral of components
               of dBperp/dxperp matrix along ray.

    switch: If true, switch the order of the integration variable and the dependent
            variable in the call(), so as to allow use with both the ode and the
            odeint interfaces.

    **kwargs: Additional parameters that may be interpreted by subclasses.

    The RHS function is implemented by the __call__(v, r) method, where

    v: ODE state vector, comprised of location, unit normal vector, J integral
       (if Jint=True), and B_perp integral (2 components) (if Bperpint=True).

    r: radial coordinate from implosion center along ray

    ---------------------------------------------
    The __call__() method will return the RHS of the ODE of the ray, integrated
    with respect to the z-direction.

    __call__() will get actual field strengths from a user-supplied function
    field(x), where x[0:2] is the location in space. field() should
    measure distances in units of l_i, and field strengths in units of b0.
    See the documentation of the stub field() below.

    The returned value is a list of floats whose length depends on Jint,
    Bperpint, and dBperpint.  return[0:3] is always the normal vector n,
    while return[3:6] is always the scaled Lorentz force.  If Jint, Bperpint,
    and dBperpint are all False, this is all that is returned.

    Jint == True, Bperpint == False, dBperpint == False:
      return is a list(7), return[6] is curl_z(B).

    Jint == False, Bperpint == True, dBperpint == False:
      return is a list(8), return[6:8] is [Bx, By].

    Jint == False, Bperpint == False, dBperpint == True:
      return is a list(9), return[6:9] is [dBx/dx, dBy/dx, dBx/dy].

    Jint == True, Bperpint == True, dBperpint == False:
      return is a list(9), return[6] is curl_z(B), return[7:9] is  [Bx, By].

    Jint == True, Bperpint == False, dBperpint == True:
      return is a list(10), return[6] is curl(B)_z, return[7:10] is [dBx/dx, dBy/dx, dBx/dy].

    Jint == False, Bperpint == True, dBperpint == True:
      return is a list(11), return[6:8] is [Bx, By], return[8:11] is [dBx/dx, dBy/dx, dBx/dy].

    Jint == True, Bperpint == True, dBperpint == True:
      return is a list(12), return[6] is curl_z(B), return[7:9] is [Bx, By], return[9:12] is [dBx/dx, dBy/dx, dBx/dy].

    Note: the z-integral of dBy/dy is -1 times the z-integral of dBx/dx, so we don't need it.
    """

#########################################
    def __init__(self, b0, li, Tkin, origin=[0,0,0], euler_angles=None, Jint=False, Bperpint=False, dBperpint=False, switch=False, **kwargs):

        self.b0 = b0
        self.li = li
        self.Tkin = Tkin
        self.origin = np.array(origin)
        self.euler_angles = euler_angles
        if euler_angles:
            self.rotmat = self.__rotation(euler_angles)

        self.Jint = Jint
        self.Bperpint = Bperpint
        self.dBperpint = dBperpint
        self.switch = switch
        self.kwargs = kwargs

        e = 4.8032E-10   # Statcoul
        mp = 1.6726E-24  # g
        c = 2.9979E+10   # cm/s
        ergperMeV = 1.6022E-06 # 1 MeV / 1 erg
        v = math.sqrt(2*(Tkin*ergperMeV)/mp)
        rlarmor = v * mp * c / (e * b0)

        self.Jconst = 1.0 / rlarmor
        self.deflangle = li / rlarmor

#########################################
    def __call__(self, v, r):

        if self.switch:
            buf = r
            r = v
            v = buf
        # Reminder: since the ODE is with respect to z, rather than with
        # respect to path length along the ray, the RHS needs to be divided
        # by n_z.

        if self.Jint and self.dBperpint:
            b, j, db = self.rotfield(v)
        elif self.Jint:
            b, j = self.rotfield(v)
        elif self.dBperpint:
            b, db = self.rotfield(v)
        else:
            b = self.rotfield(v)

        # Make sure normal vector is normalized
        n = np.array(v[3:6])
        norm = math.sqrt(np.dot(n,n))
        n = n / norm
        v[3:6] = n

        # Lorentz force
        l=[0,0,0]
        l[0] = self.deflangle * ( n[1]*b[2] - n[2]*b[1] ) / n[2]
        l[1] = self.deflangle * ( n[2]*b[0] - n[0]*b[2] ) / n[2]
        l[2] = self.deflangle * ( n[0]*b[1] - n[1]*b[0] ) / n[2]

        # Return list
        ret = list(n / n[2]) + l
        if self.Jint:
            ret.append(np.dot(n,j) / n[2])
        if self.Bperpint:
            ret = ret + list(b[0:2] / n[2])
        if self.dBperpint:
            db /= n[2]
            ret = ret + [db[0,0], db[0,1], db[1,0]]

        return ret

#########################################
    def rotfield(self, v):

        # Here we rotate the plain-vanilla field defined in field().

        # Inverse-rotate origin-shifted location
        y = v[0:3] - self.origin
        if self.euler_angles:
            x = np.dot(self.rotmat.T, y)
        else:
            x = y

        if self.Jint and self.dBperpint:
            b, j, db = self.field(x)
        elif self.Jint:
            b, j = self.field(x)
        elif self.dBperpint:
            b, db = self.field(x)
        else:
            b = self.field(x)

        # Rotate field vectors
        if self.euler_angles:
            b = np.dot(self.rotmat, b)
            if self.Jint:
                j = np.dot(self.rotmat, j)
            if self.dBperpint:
                db = np.dot(self.rotmat, db)
                db = np.dot(db, self.rotmat.T)

        if self.Jint and self.dBperpint:
            return b, j, db
        elif self.Jint:
            return b, j
        elif self.dBperpint:
            return b, db
        else:
            return b

#########################################
    def __rotation(self, euler_angles):

        c0 = math.cos(euler_angles[0])
        s0 = math.sin(euler_angles[0])
        c1 = math.cos(euler_angles[1])
        s1 = math.sin(euler_angles[1])
        c2 = math.cos(euler_angles[2])
        s2 = math.sin(euler_angles[2])

        rm = np.identity(3)
        rm[0,0] = c0
        rm[0,1] = -s0
        rm[1,0] = s0
        rm[1,1] = c0

        buf = np.identity(3)
        buf[0,0] = c1
        buf[0,2] = s1
        buf[2,0] = -s1
        buf[2,2] = c1

        rm = np.dot(buf, rm)

        buf = np.identity(3)
        buf[0,0] = c2
        buf[0,1] = -s2
        buf[1,0] = s2
        buf[1,1] = c2

        rm = np.dot(buf, rm)

        return rm

#########################################
    def field(self, x):

        # Always return b as an ndarray(3) in the first element
        # of the return tuple.
        # If self.Jint == True, add the curl(B) vector as an ndarray(3) to the
        # return tuple.
        # If self.dBperpint == True, add the matrix m[i,j]=dB_i/dx_j as an
        # ndarray((3,3)) to the return tuple.

        b = np.zeros(3)
        if self.Jint:
            j = np.zeros(3)
        if self.dBperpint:
            db = np.zeros((3,3))

        if self.Jint and self.dBperpint:
            return (b, j, db)
        elif self.Jint:
            return (b, j)
        elif self.dBperpint:
            return (b, db)
        else:
            return b

################################################################################
class pr_beercan(pr_magfield):
    """
    A subclass of pr_magfield that implements a wire-in-a-beer can-shaped
    magnetic field. The The field is azimuthally symmetric.

    The wire radius is r0, the return-current wall is at r1 and has thickness
    delta, the height of the can is z0, and the thickness of the disks at the top
    and bottom is also delta.

    B(r,z) = f(r) g(z)

            | Jr/2                              ; r < r0
            | J r0^2 /(2 r)                     ; r0 < r < r1
     f(r) = | [J r0^2 /(2 r1)] [1-(r-r1)/delta] ; r1 < r < r1+delta
            | 0                                 ; r > r1+delta

            | 1                      ; |z| < z0/2
     g(z) = | 1 - (|z|-z0/2)/delta   ; z0/2 < |z| < z0/2+delta
            | 0                      ; |z| > z0/2+delta

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['J']:  |Curl(B)| in the wire
    kwargs['r0']:  Radius of wire, in units of self.li.
    kwargs['r1']:  (Inner) radius of wall, in units of self.li.
    kwargs['delta']:  Thickness of wall and end caps.
    kwargs['z0']: Height of (interior of) can
    """

#########################################
    def field(self, x):

        # dBperpint not yet implemented
        assert not self.dBperpint

        J = self.kwargs['J']
        r0 = self.kwargs['r0']
        r1 = self.kwargs['r1']
        delta = self.kwargs['delta']
        z0 = self.kwargs['z0']

        if r0 <= 0.0 or \
           r1 <= r0 or \
           delta <= 0.0 or \
           z0 <= 0.0:
               raise ValueError

        r = math.sqrt(x[0]**2 + x[1]**2)
        cphi = x[0] / r
        sphi = x[1] / r
        z = x[2]

        if r < r0:
            fp = 0.5 * J
            f = fp*r
        elif r < r1:
            f = 0.5 * J * r0**2 / r
            fp = -f / r
        elif r < r1+delta:
            buf = 0.5 * J * r0**2 / r1
            f = buf * (1-(r-r1)/delta)
            fp = -buf/delta
        else:
            f = 0.0
            fp = 0.0

        if abs(z) < z0/2.0:
            g = 1.0
            gp = 0.0
        elif abs(z) < z0/2.0 + delta:
            g = 1 - (abs(z) - z0/2.0)/delta
            gp = - np.sign(z) / delta
        else:
            g = 0.0
            gp = 0.0

        Bphi = f * g
        b = np.zeros(3)
        b[0] = -Bphi * sphi
        b[1] =  Bphi * cphi

        if self.Jint:
            j = np.zeros(3)
            Jr = -f * gp
            Jz = g * (fp + f/r)
            j[0] = Jr * cphi
            j[1] = Jr * sphi
            j[2] = Jz
            return (b, j)
        else:
            return b

################################################################################
class pr_donut(pr_magfield):
    """
    A subclass of pr_magfield that implements a donut-shaped magnetic field. The
    axis of the donut is aligned with the z-axis.  The field is azimuthally
    symmetric.  It is tapered to zero by a poloidal current that goes linearly
    to zero in a user-set skin depth.

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['hoop_r']:  Radius of donut hoop, in units of self.li.
    kwargs['girth']:  Radius of donut girth, in units of self.li.
                        (kwargs['girth'] < kwargs['hoop_r'])
    kwargs['skin']:  Current skin depth, in units of self.li.
                        (kwargs['skin'] < kwargs['girth'])
    """

#########################################
    def field(self, x):

        # dBperpint not yet implemented
        assert not self.dBperpint

        Rh = self.kwargs['hoop_r']
        Rg = self.kwargs['girth']
        Skin = self.kwargs['skin']
        if Rh <= 0.0 or \
           Rg <= 0.0 or Rg >= Rh or \
           Skin <= 0.0 or Skin >= Rg:
               raise ValueError

        rperp = math.sqrt(x[0]**2 + x[1]**2)
        cphi = x[0] / rperp
        sphi = x[1] / rperp
        trad = math.sqrt((rperp-Rh)**2 + x[2]**2)
        irad = Rg - Skin
        Jperp = 0.0
        Jz = 0.0
        if trad >= Rg:
            Bphi = 0.0
        elif trad >= irad:
            Bphi = 1.0 - (trad - irad) / Skin
            if self.Jint:
                Jperp = x[2] / (Skin * trad)
                Jz = -(rperp - Rh) / (Skin * trad) + \
                      (1.0 - (trad - irad) / Skin) / rperp
        else:
            Bphi = 1.0

        b = np.zeros(3)
        b[0] = -Bphi * sphi
        b[1] = Bphi * cphi
        if self.Jint:
            j = np.zeros(3)
            j[0] = Jperp * cphi
            j[1] = Jperp * sphi
            j[2] = Jz
            return (b, j)
        else:
            return b

################################################################################
class pr_ellipsoidal(pr_magfield):
    """
    A subclass of pr_magfield that implements a Kugland et al.-style ellipsoidal
    blob field.

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['a']:  Radial ellipse parameter, in units of self.li.
    kwargs['b']:  Z-parameter of ellipse, in units of self.li.
    kwargs['ctr']: [x,y] list of blob center.  Defaults to [0,0]

    """

#########################################
    def field(self, x):

        a = self.kwargs["a"]
        b = self.kwargs["b"]
        if a <= 0 or b <= 0: raise ValueError
        ctr = self.kwargs.get("ctr", [0,0])

        bb = np.zeros(3)
        if self.Jint: j = np.zeros(3)
        if self.dBperpint: db = np.zeros((3,3))

        xx = x[0]-ctr[0]
        yy = x[1]-ctr[1]
        if xx != 0.0 or yy != 0.0:
            rperp = math.sqrt(xx**2 + yy**2)
            z = x[2]
            cphi = xx / rperp
            sphi = yy / rperp

            Bphi = (rperp/a) * math.exp(-(rperp/a)**2 - (z/b)**2)

            bb[0] = -Bphi * sphi
            bb[1] = Bphi * cphi
            if self.Jint:
                Jperp = 2 * z * Bphi / b**2
                Jz = 2 * Bphi * (1/rperp - rperp/a**2)
                j[0] = Jperp * cphi
                j[1] = Jperp * sphi
                j[2] = Jz

            if self.dBperpint:
                troa2 = 2 * rperp / a**2
                irp = 1.0 / rperp
                c2 = cphi**2
                s2 = sphi**2
                tzob2 = 2 * z / b**2
                db[0,0] = Bphi * cphi * sphi * troa2
                db[1,1] = -db[0,0]
                db[0,1] = Bphi * ( troa2 * s2 - irp )
                db[1,0] = Bphi * ( irp - troa2 * c2 )
                db[0,2] = Bphi * tzob2 * sphi
                db[1,2] = -Bphi * tzob2 * cphi

        elif self.Jint:
                j[2] = (2.0/a) * math.exp(-(z/b)**2)

        if self.Jint and self.dBperpint:
            return (bb, j, db)
        elif self.Jint:
            return (bb, j)
        elif self.dBperpint:
            return (bb, db)
        else:
            return bb

################################################################################
class pr_sum(pr_magfield):
    """
    A subclass of pr_magfield that performs a weighted sum of the field
    of other pr_magfield instances.

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['magfields']:  A list whose members are pr_magfield subclass
    instances, the fields and currents of which are to be summed.

    Each instance in the list is summed, weighted by its b0 parameter.
    """

#########################################
    def field(self, x):

        b = np.zeros(3)
        if self.Jint: j = np.zeros(3)
        if self.dBperpint: db = np.zeros((3,3))

        for magfield in self.kwargs['magfields']:
            magfield.Jint = self.Jint
            magfield.dBperpint = self.dBperpint
            buf = magfield.rotfield(x)
            b += buf[0] * magfield.b0
            if self.Jint and self.dBperpint:
                j += buf[1] * magfield.b0
                db += buf[2] * magfield.b0
            elif self.Jint:
                j += buf[1] * magfield.b0
            elif self.dBperpint:
                db += buf[1] * magfield.b0

        if self.Jint and self.dBperpint:
            return (b, j, db)
        elif self.Jint:
            return (b,j)
        elif self.dBperpint:
            return (b, db)
        else:
            return b

################################################################################
class pr_mesh(pr_magfield):
    """
    A subclass of pr_magfield that reads in a cubical 3-dimensional mesh of
    magnetic field values from a file at initialization, and returns values
    from that mesh.

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['meshfile']:  A filename to be opened.  The file should contain the
    mesh values.  The logical size of the cube will be read in from the file's
    header data.  If required, the mesh values of the 3-d current vector and
    3x3 matrix of field derivatives will also be expected to be available in
    the file. See the constructor for details.

    The physical dimension of the sides of the mesh is self.li, so of course
    in the units adopted in this class family the numeric value of the side length
    is 1.  The origin is at the center of the mesh.

    The default file format is a very basic ASCII table with a simple header.
    To read more complex formats, subclass this class and override the readmesh()
    method.
    """

#########################################
    # Overriding constructor to read in file at initialization.
    #
    def __init__(self, b0, li, Tkin, origin=[0,0,0], euler_angles=None, Jint=False, Bperpint=False, dBperpint=False, switch=False, **kwargs):

        super(pr_mesh, self).__init__(b0, li, Tkin, origin, euler_angles, Jint, Bperpint, dBperpint, switch, **kwargs)

        mf = open(kwargs["meshfile"], 'r')

        self.readmesh(mf)

#########################################
    def readmesh(self, mf):

#        print "Reading Mesh File..."
        line = mf.readline()
        while re.search('^#', line):

            if re.search('=', line):
                key = line[1:].split("=")[0].strip()
                val = line[1:].split("=")[1].strip()
                self.kwargs[key] = val

            if re.search('^# Cells per dimension = ', line):
                self.Ncell = int(line.split()[5])
                self.delta = 1.0 / self.Ncell
            line = mf.readline()

        self.Bmesh = np.zeros((self.Ncell, self.Ncell, self.Ncell, 3))
        dbp0 = 6
        if self.Jint:
            self.Jmesh = np.zeros((self.Ncell, self.Ncell, self.Ncell, 3))
            dbp0 += 3
        if self.dBperpint:
            self.dbmesh = np.zeros((self.Ncell, self.Ncell, self.Ncell, 3, 3))

        while line:
            ix = int(line.split()[0])
            iy = int(line.split()[1])
            iz = int(line.split()[2])
            self.Bmesh[ix, iy, iz, 0] = float(line.split()[3])
            self.Bmesh[ix, iy, iz, 1] = float(line.split()[4])
            self.Bmesh[ix, iy, iz, 2] = float(line.split()[5])
            if self.Jint:
                self.Jmesh[ix, iy, iz, 0] = float(line.split()[6])
                self.Jmesh[ix, iy, iz, 1] = float(line.split()[7])
                self.Jmesh[ix, iy, iz, 2] = float(line.split()[8])
            if self.dBperpint:
                self.dbmesh[ix, iy, iz, 0, 0] = float(line.split()[dbp0])
                self.dbmesh[ix, iy, iz, 0, 1] = float(line.split()[dbp0+1])
                self.dbmesh[ix, iy, iz, 0, 2] = float(line.split()[dbp0+2])
                self.dbmesh[ix, iy, iz, 1, 0] = float(line.split()[dbp0+3])
                self.dbmesh[ix, iy, iz, 1, 1] = float(line.split()[dbp0+4])
                self.dbmesh[ix, iy, iz, 1, 2] = float(line.split()[dbp0+5])
                self.dbmesh[ix, iy, iz, 2, 0] = float(line.split()[dbp0+6])
                self.dbmesh[ix, iy, iz, 2, 1] = float(line.split()[dbp0+7])
                self.dbmesh[ix, iy, iz, 2, 2] = float(line.split()[dbp0+8])

            line = mf.readline()

        mf.close()
#        print "...done"

#########################################
    def field(self, x):

        b = np.zeros(3)
        if self.Jint: j = np.zeros(3)
        if self.dBperpint: db = np.zeros((3,3))

        if x[0] >= -0.5 and x[0] <= 0.5 and x[1] >= -0.5 and x[1] <= 0.5 and x[2] >= -0.5 and x[2] <= 0.5:
            ix = min(int((x[0] + 0.5) / self.delta), self.Ncell-1)
            iy = min(int((x[1] + 0.5) / self.delta), self.Ncell-1)
            iz = min(int((x[2] + 0.5) / self.delta), self.Ncell-1)
            b = self.Bmesh[ix,iy,iz,:]
            if self.Jint: j = self.Jmesh[ix,iy,iz,:]
            if self.dBperpint: db = self.dbmesh[ix,iy,iz,:,:]

        if self.Jint and self.dBperpint:
            return (b, j, db)
        elif self.Jint:
            return (b,j)
        elif self.dBperpint:
            return (b, db)
        else:
            return b


################################################################################
class pr_mesh_gp(pr_mesh):
    """
    A subclass of pr_mesh that reads in a cubical 3-dimensional mesh of
    Gaussian process interpolant weights from a file, and returns interpolated
    field values.

    The self.kwargs dictionary read in by the superclass constructor is interpreted as
    follows:

    kwargs['meshfile']:  A filename to be opened.  The file should contain the
    mesh weight values.  The logical size of the cube will be read in from the file's
    header data. See the constructor for details.

    The physical dimension of the sides of the mesh is self.li, so of course
    in the units adopted in this class family the numeric value of the side length
    is 1.  The origin is at the center of the mesh.

    The default file format is a very basic ASCII table with a simple header.
    To read more complex formats, subclass this class and override the readmesh()
    method.
    """

#########################################
    def readmesh(self, mf):

#        print "Reading Mesh File..."
        stcl = []
        line = mf.readline()
        while re.search('^#', line):

            if re.search('=', line):
                key = line[1:].split("=")[0].strip()
                val = line[1:].split("=")[1].strip()
                self.kwargs[key] = val

            if re.search('# SE_sigma = ', line):
                sigma = float(line.split()[3])
                self.sm2 = 1.0 / sigma**2
                self.sm4 = self.sm2**2
                self.sm6 = self.sm2**3

            if re.search('# GP Stencil point: ', line):
                p = [line.split()[4], line.split()[5], line.split()[6]]
                stcl.append(p)

            if re.search('^# Cells per dimension = ', line):
                self.Ncell = int(line.split()[5])
                self.delta = 1.0 / self.Ncell

            line = mf.readline()

        self.stencil = np.array(stcl, dtype=np.int_)
        Nstencil = len(self.stencil)
        self.Warr = np.zeros((self.Ncell, self.Ncell, self.Ncell, 3*Nstencil))

        while line:
            l = line.split()
            ix = int(l[0])
            iy = int(l[1])
            iz = int(l[2])
            self.Warr[ix, iy, iz, :] = l[9:3*Nstencil+9]
            line = mf.readline()

        mf.close()

#########################################
    def field(self, x):

        # dBperpint not yet implemented
        assert not self.dBperpint

        Nstencil = len(self.stencil)
        b = np.zeros(3)
        bv = np.zeros((3, Nstencil, 3))
        if self.Jint:
            j = np.zeros(3)
            jv = np.zeros((3, Nstencil, 3))

        if x[0] >= -0.5 and x[0] <= 0.5 and x[1] >= -0.5 and x[1] <= 0.5 and x[2] >= -0.5 and x[2] <= 0.5:
            xx = (x[0] + 0.5) / self.delta #
            yy = (x[1] + 0.5) / self.delta # Note unit change from li to delta!
            zz = (x[2] + 0.5) / self.delta #
            ix = min( int(xx), self.Ncell-1 )
            iy = min( int(yy), self.Ncell-1 )
            iz = min( int(zz), self.Ncell-1 )
            for ist in range(Nstencil):
                pp = self.stencil[ist] + [ix+0.5, iy+0.5, iz+0.5]
                dx = np.array([xx, yy, zz]) - pp
                dx2 = dx.dot(dx)
                e = math.exp(-0.5 * dx2 * self.sm2)
                sm2e = self.sm2 * e
                sm4e = self.sm4 * e
                for l in range(3):
                    bv[l, ist, l] = 2*sm2e - dx2*sm4e
                    for lp in range(3):
                        bv[l, ist, lp] += sm4e * dx[l] * dx[lp]

                if self.Jint:
                    jf = (5*self.sm4 - dx2*self.sm6) * e / self.delta # 1/delta is consequence of unit change
                    jv[0, ist, 1] = jf * dx[2]
                    jv[1, ist, 2] = jf * dx[0]
                    jv[2, ist, 0] = jf * dx[1]
                    jv[1, ist, 0] = -jv[0, ist, 1]
                    jv[2, ist, 1] = -jv[1, ist, 2]
                    jv[0, ist, 2] = -jv[2, ist, 0]

            bvv = bv.reshape((3, 3*Nstencil))
            b = bvv.dot(self.Warr[ix,iy,iz,:])
            if self.Jint:
                jvv = jv.reshape((3, 3*Nstencil))
                j = jvv.dot(self.Warr[ix,iy,iz,:])

        if self.Jint:
            return(b,j)
        else:
            return b


################################################################################
class pr_trace(object):

    """
    Ray integration class.

    The idea is that the magnetized domain mag -- an instance of a subclass of
    pr_magfield -- is placed with its origin at the coordinates x=0, y=0, z=ri
    with respect to the experimental coordinates.  In experimental coordinates,
    the source is at the origin and the detector is an infinite plane normal
    to the z-axis, located at z=rs.

    A circular aperture of radius raperture, with a plane normal to the
    z-direction, and with its center at x=0, y=0, z=ri, is the target for all
    proton rays.

    The __call__ invocation is used to request a number of randomly-generated
    rays.

    -------------------------------
    Constructor invocation:

    __init__(self, mag, rs, ri, raperture, rng_seed=None, tol=1.49012e-8, h0=0.1, hmax=0.1, Nstep=1)

    mag:       An instance of a subclass of pr_magfield, specifying the magnetized
               domain

    rs:        Distance along the z-axist to the detector screen from the origin

    ri:        z-coordinate of the center of the magnetized domain

    raperture: Radius of target aperture

    rng_seed:  Optional, resets numpy's RNG seed.

    tol:       ODE tolerance -- see documentation for rtol and atol in
               scipy.integrate.odeint

    h0:        Initial step for ODE solver -- see documentation for h0 in
               scipy.integrate.odeint

    hmax:      Largest allowed step for ODE solver -- see documentation for
               hmax in scipy.integrate.odeint

    hmin:      Smallest allowed step for ODE solver -- see documentation for
               hmax in scipy.integrate.odeint

    Nstep:     Number of steps to take across domain.

    mxordn:    Max order for non-stiff integrator -- see documentation for
               hmax in scipy.integrate.odeint

    -------------------------------
    Call invocation:

    __call__(self, nray)

    nray:    Number of rays to be simulated

    The output is a list, each element of which comprises the data for a single
    ray.  The data for a ray include
      -- the initial normal vector (a 3-D list)
      -- the location of the final position of the ray on the screen (a 2-D list)
      -- the integrated current component along the ray, if mag.Jint=True
      -- the integrated perpendicular magnetic field (a 2-D list (Bx, By)), if
          mag.Bperpint=True

    -------------------------------
    Other methods available:

    -------------------------------
    print_header(self, fd)

    Prints a commented header containing the parameters of the experiment to
    file descriptor fd

    fd: File descriptor, opened by user.

    -------------------------------
    print_data(self, fd, data)

    Prints the data, output by the __call__(), to the file descriptor fd.

    fd: File descriptor, opened by user.
    data: data array, output by __call__()

    """

#########################################
    def __init__(self, mag, rs, ri, raperture, rng_seed=32344541, tol=1.49012e-8, h0=0.1, hmax=0.1, hmin=0.0, Nstep=1, mxordn=12):

        self.mag = mag
        self.rs = rs
        self.ri = ri
        self.raperture = raperture
        self.tol = tol
        self.rng_seed = rng_seed
        self.h0 = h0
        self.hmax = hmax
        self.hmin = hmin
        self.Nstep = Nstep
        self.mxordn = mxordn
        np.random.seed(rng_seed)

#########################################
    def __call__(self, nray):

        results = []
        for iray in range(nray):

            line = []
            nrm = self.initnormal() # Initial normal vector of ray
            line += nrm
            n2 = nrm[2]

            z = self.ri - self.mag.li / 2.0
            x = z * (nrm[0] / nrm[2]) / self.mag.li
            y = z * (nrm[1] / nrm[2]) / self.mag.li

            initstate = [x, y, -0.5] + nrm
            bp0 = 6
            dbp0 = 6
            if self.mag.Jint:
                initstate += [0.0]
                bp0 += 1
                dbp0 += 1
            if self.mag.Bperpint:
                initstate += [0.0,0.0]
                dbp0 += 2
            if self.mag.dBperpint:
                initstate += [0.0,0.0,0.0]

            # Integrate ray in domain, from z=-li/2 to z=+li/2, which is
            # z=-0.5 to z=+0.5 in the coordinates used in the domain.
            zpts = np.linspace(-0.5, 0.5, num=self.Nstep+1)
            exitpt= odeint(self.mag, initstate, zpts, rtol=self.tol, \
                            atol=self.tol, h0 = self.h0, hmax = self.hmax, hmin = self.hmin, mxordn=self.mxordn)

            #exitpt, fo = odeint(self.mag, initstate, zpts, rtol=self.tol, \
                            #atol=self.tol, h0 = self.h0, hmax = self.hmax, hmin = self.hmin, \
                            #full_output=True)
            #foo=open('foo','w')
            #for it in fo.items(): foo.write("%s : %s\n" % it)
            #foo.close

            z = self.ri + self.mag.li / 2.0
            dz = self.rs - z
            x = exitpt[self.Nstep][0] * self.mag.li
            y = exitpt[self.Nstep][1] * self.mag.li
            npnrm = exitpt[self.Nstep][3:6]
            nrm = list(npnrm / math.sqrt(npnrm.dot(npnrm)))
            xf = x + dz * nrm[0] / nrm[2]
            yf = y + dz * nrm[1] / nrm[2]

            line += [xf, yf]
            line += nrm
            if self.mag.Jint:
                xi = exitpt[self.Nstep][6] * self.mag.Jconst * self.ri \
                                    * (self.rs - self.ri) / self.rs
                xi *= n2 # d/ds, not d/dz
                line.append(xi)

            if self.mag.Bperpint: line += list(exitpt[self.Nstep][bp0:bp0+2] * (self.mag.b0*self.mag.li*n2))

            if self.mag.dBperpint: line += list(exitpt[self.Nstep][dbp0:dbp0+3] * self.mag.b0 * n2)

            results.append(line)

        return results

#########################################
    def initnormal(self):

        c0 = self.ri / math.sqrt(self.raperture**2 + self.ri**2)
        rx = uniform.rvs(size=2)
        ctheta = rx[0] * (1.0 - c0) + c0
        stheta = math.sqrt(1.0 - ctheta**2)
        phi = rx[1] * 2.0 * math.pi

        nx = stheta * math.cos(phi)
        ny = stheta * math.sin(phi)
        nz = ctheta

        return [nx, ny, nz]

#########################################
    def print_header(self, fd):

        fd.write('# Number of MPI ranks: %d\n' % MPI.COMM_WORLD.size)
        fd.write('# Field Configuration Type: ' + type(self.mag).__name__ + '\n')
        fd.write('# Field Configuration Parameters: ' + str(self.mag.kwargs) + '\n')
        fd.write('# B0: %12.5E (Gauss)\n' % self.mag.b0)
        fd.write('# li: %12.5E (cm)\n' % self.mag.li)
        fd.write('# Tkin: %12.5E (MeV)\n' % self.mag.Tkin)
        fd.write('# Euler Angles: ' + str(self.mag.euler_angles) + '\n')
        fd.write('# Jint: ' + str(self.mag.Jint) + '\n')
        fd.write('# Bperpint: ' + str(self.mag.Bperpint) + '\n')
        fd.write('#\n# rs: %12.5E (cm)\n' % self.rs)
        fd.write('# ri: %12.5E (cm)\n' % self.ri)
        fd.write('# raperture: %12.5E (cm)\n' % self.raperture)
        fd.write('# rng_seed: %d\n' % self.rng_seed)
        fd.write('# tol: %12.5E\n' % self.tol)
        fd.write('# Initial step size (h0): %12.5E\n' % self.h0)
        fd.write('# Max step size (hmax): %12.5E\n' % self.hmax)
        fd.write('#\n# Columns:\n#\n')
        fd.write('# 1: Initial normal vector X-component\n')
        fd.write('# 2: Initial normal vector Y-component\n')
        fd.write('# 3: Initial normal vector Z-component\n')
        fd.write('# 4: Final X location at screen (cm)\n')
        fd.write('# 5: Final Y location at screen (cm)\n')
        fd.write('# 6: Final normal vector X-component\n')
        fd.write('# 7: Final normal vector Y-component\n')
        fd.write('# 8: Final normal vector Z-component\n')
        nb0 = 9
        if self.mag.Jint:
            fd.write('# 9: Current integral\n')
            nb0 += 1
        if self.mag.Bperpint:
            fd.write('# %1d: Bx Integral\n' % nb0)
            fd.write('# %1d: By Integral\n' % (nb0+1))
            nb0 += 2
        if self.mag.dBperpint:
            fd.write('# %1d: dBx/dx Integral\n' % nb0)
            fd.write('# %1d: dBx/dy Integral\n' % (nb0+1))
            fd.write('# %1d: dBy/dx Integral\n' % (nb0+2))
        fd.write('#\n')

#########################################
    def print_data(self, fd, data):

        for line in data:
            s = ''
            for datum in line: s += '%12.5E  ' % datum
            s += '\n'
            fd.write(s)


################################################################################
class pr_trace_mpi(pr_trace):

    """
    Ray integration class, MPI-aware. Construction and invocation are the
    same as for the pr_trace parent class.  The difference is that __init__()
    initializes the RNG using self.random_seed+rank, where rank is the MPI
    rank, and __call__() divides up the number of protons between the available
    ranks.
    """

#########################################
    def __init__(self, mag, rs, ri, raperture, rng_seed=32344541, tol=1.49012e-8, h0=0.1, hmax=0.1, hmin=0.0, Nstep=1, mxordn=12):

        super(pr_trace_mpi, self).__init__(mag, rs, ri, raperture, rng_seed, tol, h0, hmax, hmin, Nstep, mxordn)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        np.random.seed(rng_seed + rank)


#########################################
    def __call__(self, nray):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        rays_per_rank = nray / size
        remainder = nray % size

        if rank == 0:
            results = super(pr_trace_mpi, self).__call__(rays_per_rank + remainder)
            for ir in range(1, size):
                results += comm.recv(source=ir, tag=10+ir)
            return results
        else:
            results = super(pr_trace_mpi, self).__call__(rays_per_rank)
            comm.send(results, dest=0, tag=10+rank)

################################################################################
class pr_trace2(pr_trace):

    """
    New, improved ray integration class, simplified interface, and with MPI baked in.

    As before, the magnetized domain mag -- an instance of a subclass of
    pr_magfield -- is placed with its origin at the coordinates x=0, y=0, z=ri
    with respect to the experimental coordinates.  In experimental coordinates,
    the source is at the origin and the detector is an infinite plane normal
    to the z-axis, located at z=rs.

    A circular aperture of radius raperture, with a plane normal to the
    z-direction, and with its center at x=0, y=0, z=ri, is the target for all
    proton rays.

    The __call__ invocation is used to request a number of randomly-generated
    rays.

    The interface simplification occurs because all ODE parameters are replaced
    by an instance of ode from scipy.integrate, which is prepared by the calling
    program.  The magnetic field instance of a subclass of pr_magfield() should
    be embedded in the ode instance, as the RHS of the ODE.

    -------------------------------
    Constructor invocation:

    __init__(self, rs, ri, raperture, integrator, Nstep=1, rng_seed=None)

    rs:        Distance along the z-axist to the detector screen from the origin

    ri:        z-coordinate of the center of the magnetized domain

    raperture: Radius of target aperture

    odeinst: Instance of ode from scipy.integrate, set from calling program.

    Nstep: Number of steps to take in integration.

    rng_seed:  Optional, resets numpy's RNG seed. Will be summed to the MPI rank
               of each mpi process.

    """

#########################################
    def __init__(self, rs, ri, raperture, odeinst, Nstep=1, rng_seed=None):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        self.rng_seed = rng_seed
        np.random.seed(rng_seed + rank)
        self.rs = rs
        self.ri = ri
        self.raperture = raperture
        self.odeinst = odeinst
        self.Nstep = Nstep
        self.mag = odeinst.__dict__['f']
        self.tol = odeinst._integrator.__dict__['rtol']
        self.h0 = odeinst._integrator.__dict__['first_step']
        self.hmax = odeinst._integrator.__dict__['max_step']

#########################################
    def __call__(self, nray):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        rays_per_rank = nray / size
        remainder = nray % size

        if rank == 0:
            results = self.__integrate_ray(rays_per_rank + remainder)
            for ir in range(1, size):
                results += comm.recv(source=ir, tag=10+ir)
            return results
        else:
            results = self.__integrate_ray(rays_per_rank)
            comm.send(results, dest=0, tag=10+rank)


 #########################################
    def __integrate_ray(self, nray):

        results = []
        for iray in range(nray):

            line = []
            nrm = self.initnormal() # Initial normal vector of ray
            line += nrm
            n2 = nrm[2]

            z = self.ri - self.mag.li / 2.0
            x = z * (nrm[0] / nrm[2]) / self.mag.li
            y = z * (nrm[1] / nrm[2]) / self.mag.li

            initstate = [x, y, -0.5] + nrm
            bp0 = 6
            dbp0 = 6
            if self.mag.Jint:
                initstate += [0.0]
                bp0 += 1
                dbp0 += 1
            if self.mag.Bperpint:
                initstate += [0.0,0.0]
                dbp0 += 2
            if self.mag.dBperpint:
                initstate += [0.0,0.0,0.0]

            # Integrate ray in domain, from z=-li/2 to z=+li/2, which is
            # z=-0.5 to z=+0.5 in the coordinates used in the domain.
            self.odeinst.set_initial_value(initstate, -0.5)
            exitpt = np.zeros((self.Nstep+1, len(initstate)))
            exitpt[0] = initstate
            zpts = np.linspace(-0.5, 0.5, num=self.Nstep+1)
            for step in range(1, self.Nstep+1):
                exitpt[step] = self.odeinst.integrate(zpts[step])
                if not self.odeinst.successful(): sys.exit("Integration Failure")

            z = self.ri + self.mag.li / 2.0
            dz = self.rs - z
            x = exitpt[self.Nstep][0] * self.mag.li
            y = exitpt[self.Nstep][1] * self.mag.li
            npnrm = exitpt[self.Nstep][3:6]
            nrm = list(npnrm / math.sqrt(npnrm.dot(npnrm)))
            xf = x + dz * nrm[0] / nrm[2]
            yf = y + dz * nrm[1] / nrm[2]

            line += [xf, yf]
            line += nrm
            if self.mag.Jint:
                xi = exitpt[self.Nstep][6] * self.mag.Jconst * self.ri \
                                    * (self.rs - self.ri) / self.rs
                xi *= n2 # d/ds, not d/dz
                line.append(xi)

            if self.mag.Bperpint: line += list(exitpt[self.Nstep][bp0:bp0+2] * (self.mag.b0*self.mag.li*n2))

            if self.mag.dBperpint: line += list(exitpt[self.Nstep][dbp0:dbp0+3] * self.mag.b0 * n2)

            results.append(line)

        return results
