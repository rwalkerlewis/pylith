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
# @file tests/fullscale/linearporoelasticity/nofaults-2d/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
#  
#        -1000 Pa
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         Uy=0
#
# Dirichlet boundary conditions
#   Ux(+-4000,0) = 0
#   Uy(x,-4000) = 0
# Neumann boundary conditions
#   \tau_normal(x,+4000) = -1000
#

import numpy


# Physical properties
p_porosity               =      0.1 # frac
p_density                =   2500.0 # kg/m**3
p_fluid_density          =     1000 # kg/m**3
p_fluid_viscosity        =    0.001 # Pa*s
p_bulk_modulus           =    0.7e9 # Pa
p_shear_modulus          =    0.4e9 # Pa
p_biot_coefficient       =      0.8 # -
p_isotropic_permeability = 10.0e-14 # m**2
p_fluid_bulk_modulus     =      2e9 # Pa
p_external_force         =     1000 # Pa

p_bulk_density = (1.0 - p_porosity)*p_density + p_porosity*p_fluid_density
p_drained_bulk_modulus = -1.0*p_bulk_modulus * (p_biot_coefficient - 1.0)


p_mu = p_shear_modulus
p_lambda = p_drained_bulk_modulus - 2.0/3.0 * p_shear_modulus

gacc = 9.80665  # m/s
ymax = +4000.0
ymin = -4000.0  # m

# Height of column, m
H = ymax - ymin

# Drained Matrix Bulk Modulus, Pa
K_dr = (-1*p_bulk_modulus)*(p_biot_coefficient - 1)

# Confined compressibility of the Porous Medium
m_v = 1 / (K_dr + 4/3 * p_shear_modulus)

# Drained storage coefficient
S = (p_biot_coefficient - p_porosity)/K_sg + p_porosity/p_fluid_bulk_modulus

# Consolidation coefficient
c_v = p_isotropic_permeability / (p_fluid_viscosity * (p_biot_coefficient**2 * m_v + S))

# Undrained Fluid Pressure Response
p_0 = (p_biot_coefficient*m_v)/(p_biot_coefficient**2*m_v+S) * p_external_force

# Undrained response of consolidation
u_z0 = (m_v - (p_biot_coefficient**2*m_v**2)/(p_biot_coefficient**2*m_v + S))*p_external_force*H

# Steady State Consolidation
u_z_inf = m_v*p_external_force*H 


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to gravitational body forces (no initial stress).
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "gravitational_acceleration": self.gacc,
            "initial_amplitude": self.zero_vector,
        }
        return

    def getField(self, name, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        strain = self.strain(locs)

        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[:, :, 1] = p_bulk_density * gacc / (p_lambda + 2 * p_mu) * \
            (0.5 * (locs[:, 1]**2 - ymin**2) - ymax * (locs[:, 1] - ymin))
        return disp

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def density(self, locs):
        """
        Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_mu * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def bulk_modulus(self, locs):
        """
        Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = (p_lambda + 2.0 / 3.0 * p_mu) * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        eyy = p_density * gacc * (locs[:, 1] - ymax) / (p_lambda + 2 * p_mu)
        exx = 0
        ezz = 0
        exy = 0

        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = exx
        strain[:, :, 1] = eyy
        strain[:, :, 2] = ezz
        strain[:, :, 3] = exy
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        syy = p_density * gacc * (locs[:, 1] - ymax)
        sxx = p_lambda / (p_lambda + 2 * p_mu) * syy
        sxy = 0.0
        szz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sxx + syy)

        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = sxy
        return stress

    def gacc(self, locs):
        """Compute gravitational acceleration at locations.
        """
        (npts, dim) = locs.shape
        gravacc = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        gravacc[0, :, 1] = -gacc
        return gravacc


# End of file
