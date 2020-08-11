# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/linearporoelasticity/terzaghi/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
#
#          1 Pa
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         Uy=0
#
# Dirichlet boundary conditions
#   Ux(+-1,0) = 0
#   Uy(x,-1) = 0
# Neumann boundary conditions
#   \tau_normal(x,+1) = 1*Pa

import numpy

# Physical properties
p_shear_modulus = 3.0 # Pa
p_undrained_bulk_modulus = 9.76 # Pa
p_biot_coefficient = 0.6 # -
p_biot_modulus = 16.0 # Pa
p_isotropic_permeability = 1.5 # m**2
p_fluid_dynamic_viscosity = 1.0 # Pa*s
p_zmax = 1.0 # m
p_zmin = 0.0 # m
p_ymax = 10.0 # m
p_ymin = 0.0 # m
p_xmax = 10.0 # m
p_xmin = 0.0 # m
p_vertical_stress = 1.0 # Pa

# Height of column, m
H = zmax - zmin
L = H

p_drained_bulk_modulus = p_undrained_bulk_modulus - p_biot_coefficient*p_biot_coefficient*p_biot_modulus # Pa,      Cheng (B.5)
p_nu = (3.0*p_drained_bulk_modulus - 2.0*p_shear_modulus) / (2.0*(3.0*p_drained_bulk_modulus + p_shear_modulus)) # -,       Cheng (B.8)
p_nu_u = (3.0*p_undrained_bulk_modulus - 2.0*p_shear_modulus) / (2.0*(3.0*p_undrained_bulk_modulus + p_shear_modulus)) # -,       Cheng (B.9)
p_eta = (3.0*p_biot_coefficient*p_shear_modulus) /(3.0*p_drained_bulk_modulus + 4.0*p_shear_modulus) #  -,       Cheng (B.11)
p_S = (3.0*p_undrained_bulk_modulus + 4.0*p_shear_modulus) / (p_biot_modulus*(3.0*p_drained_bulk_modulus + 4.0*p_shear_modulus)) # Pa^{-1}, Cheng (B.14)
p_c = (p_isotropic_permeability / p_fluid_dynamic_viscosity) / p_S # m^2 / s, Cheng (B.16)

# Time steps
tsteps = numpy.arange(1.0, 5.01, 1.0) + 1 # sec

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to Terzaghi's problem
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 200

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        z = locs[:,1] - zmax
        t_track = 0

        for t in tsteps:
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            displacement[t_track,:,1] = ( (p_vertical_stress*L*(1.0 - 2.0*p_nu_u)) / (2.0*p_shear_modulus*(1.0 - p_nu_u)) ) \
                                        * (1.0 - z_star) + ((p_vertical_stress*L*(p_nu_u - p_nu)) / (2.0*p_shear_modulus*(1.0 - p_nu_u)*(1.0 - p_nu)))*self.F2(z_star, t_star)
            t_track += 1

        return displacement

    def pressure(self, locs):
        """
        Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - zmax
        t_track = 0

        for t in tsteps:
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            pressure[t_track,:] = (p_vertical_stress * p_eta) / (p_shear_modulus * p_S) * self.F1(z_star, t_star)
            t_track += 1

        return pressure

    def trace_strain(self, locs):
        """
        Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - zmax
        t_track = 0

        for t in tsteps:
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            trace_strain[t_track,:] = -((p_vertical_stress*L*(1.0 - 2.0*p_nu_u)) / (2.0*p_shear_modulus*(1.0 - p_nu_u)*L)) \
                                      + ((p_vertical_stress*L*(p_nu_u - p_nu)) / (2.0*p_shear_modulus*(1.0 - p_nu_u)*(1.0 - p_nu)))*self.F3(z_star, t_star)
            t_track += 1

        return trace_strain


    # Series functions

    def F1(self, z_star, t_star):
        F1 = 0
        for m in numpy.arange(1,self.ITERATIONS*2+1,2):
            F1 += 4/(m*numpy.pi) * numpy.sin(0.5*m*numpy.pi*z_star)*numpy.exp(-(m*numpy.pi)**2*t_star)
        return F1

    def F2(self, z_star, t_star):
        F2 = 0
        for m in numpy.arange(1,self.ITERATIONS*2+1,2):
            F2 += 8/(m*numpy.pi)**2  * numpy.cos(0.5*m*numpy.pi*z_star) * (1 - numpy.exp(-(m * numpy.pi)**2 * t_star) )
        return F2

    def F3(self, z_star, t_star):
        F3 = 0
        for m in numpy.arange(1, self.ITERATIONS*2+1,2):
            F3 += (-4.0 / (m*numpy.pi*L)) * numpy.sin(0.5*m*numpy.pi*z_star) * (1.0 - numpy.exp(-(m*numpy.pi)*t_star))
        return F3

# End of file
