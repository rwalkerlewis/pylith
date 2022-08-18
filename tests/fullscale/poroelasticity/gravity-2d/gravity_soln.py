# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/poroelasticity/gravity-2d/gravity_soln.py
#
# @brief Analytical solution to gravitational body forces (no initial stress).
# 
#          P=0
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         Uy=0
#
# Dirichlet boundary conditions
# Ux(+-4000,0) = 0
# Uy(x,-4000) = 0

import numpy


# Physical properties
p_solid_density = 2500.0  # kg/m**3
p_fluid_density = 1000.0 # kg/m**3
p_fluid_viscosity = 0.001 # Pa*s
p_porosity = 0.1 # -
p_shear_modulus = 6e9 # Pa
p_drained_bulk_modulus = 10e9 # Pa
p_fluid_bulk_modulus = 2e9 # Pa
p_solid_bulk_modulus = 20e9 # Pa
p_biot_coefficient = 1.0 # -
p_isotropic_permeability = 1e-14 # m**2

p_biot_modulus = 1.0 / ( p_porosity / p_fluid_bulk_modulus + (p_biot_coefficient - p_porosity) / p_solid_bulk_modulus )
p_bulk_density = p_solid_density*(1.0 - p_porosity) + p_fluid_density * p_porosity

p_mu = p_shear_modulus
p_undrained_bulk_modulus = p_drained_bulk_modulus + p_biot_coefficient**2 + p_biot_modulus
p_drained_lambda = p_drained_bulk_modulus - (2.0 * p_shear_modulus) / 3.0
p_undrained_lambda = p_undrained_bulk_modulus - (2.0 * p_shear_modulus) / 3.0
p_skempton_coefficient = (p_biot_coefficient * p_biot_modulus) / p_undrained_bulk_modulus
p_constant_stress_storage = p_undrained_bulk_modulus / (p_biot_modulus * (p_undrained_bulk_modulus - p_biot_coefficient**2 * p_biot_modulus))
p_drained_poisson_ratio = p_drained_lambda / (2.0 * (p_drained_lambda + p_mu))
p_undrained_poisson_ratio = p_undrained_lambda / (2.0 * (p_undrained_lambda + p_mu))
p_drained_youngs_modulus = 2.0 * p_shear_modulus * (1.0 * p_drained_poisson_ratio)
p_undrained_youngs_modulus = 2.0 * p_shear_modulus * (1.0 * p_undrained_poisson_ratio)

gacc = 9.80665  # m/s
ymax = +4000.0
ymin = -4000.0  # m

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to gravitational body forces (no initial stress).
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.zero_scalar,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,            
            "porosity": self.porosity,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": self.zero_vector,
        }
        return

    def getField(self, name, mesh_entity, pts):
        return self.fields[name](pts)

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        strain = self.strain(locs)

        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        # disp[:,:, 1] = ((1 + p_drained_poisson_ratio)*(1-2.0*p_drained_poisson_ratio))/(2.0*p_drained_youngs_modulus*(1.0 - p_drained_poisson_ratio)) * (-p_biot_coefficient*p_fluid_density + p_bulk_density) * gacc * (locs[:, 1] - ymin) * ((locs[:, 1] - ymin) - 2.0*(ymax - ymin))
        disp[:,:, 1] = p_bulk_density * gacc / (p_undrained_lambda + 2 * p_mu) * \
            (0.5 * (locs[:, 1]**2 - ymin**2) - ymax * (locs[:, 1] - ymin))
        return disp

    def pressure(self, locs):
        """ Assume hydrostatic pressure.
        """
        (npts, dim) = locs.shape
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        (npts, dim) = locs.shape
        pressure[0,:,0] = -p_fluid_density * gacc * (locs[:, 1] - ymax)
        return pressure

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, 1), dtype=numpy.float64)

    def solid_density(self, locs):
        """Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        solid_density = p_solid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_shear_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def porosity(self, locs):
        """Compute porosity field at locations.
        """
        (npts, dim) = locs.shape
        porosity = p_porosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def fluid_viscosity(self, locs):
        """Compute fluid_viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = p_fluid_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def drained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        undrained_bulk_modulus = p_drained_bulk_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

    def biot_coefficient(self, locs):
        """Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = p_biot_coefficient * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """Compute biot modulus field at locations.
        """
        (npts, dim) = locs.shape
        biot_modulus = p_biot_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = p_isotropic_permeability * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def strain(self, locs):
        """Compute strain field at locations.
        """
        # pressure = -p_fluid_density * gacc * (locs[:, 1] - ymax)
        # syy = p_bulk_density * gacc * (locs[:, 1] - ymax)
        # eyy = (p_biot_coefficient * pressure + syy) / (p_drained_lambda + 2.0 * p_shear_modulus)
        # exx = 0
        # ezz = 0
        # exy = 0

        eyy = p_bulk_density * gacc * (locs[:, 1] - ymax) / (p_undrained_lambda + 2 * p_mu)
        exx = 0
        ezz = 0
        exy = 0

        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:,:, 0] = exx
        strain[:,:, 1] = eyy
        strain[:,:, 2] = ezz
        strain[:,:, 3] = exy
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        pressure = -p_fluid_density * gacc * (locs[:, 1] - ymax)

        # syy = p_bulk_density * gacc * (locs[:, 1] - ymax)
        # eyy = (p_biot_coefficient * pressure + syy) / (p_drained_lambda + 2.0 * p_shear_modulus)
        # sxx = p_drained_lambda * eyy - p_biot_coefficient * pressure
        # # szz = p_drained_poisson_ratio * (sxx + syy + 2.0 * p_biot_coefficient * pressure) - p_biot_coefficient * pressure
        # szz = p_drained_poisson_ratio * (sxx + syy) - p_biot_coefficient * (1.0 - 2.0 * p_drained_poisson_ratio) * pressure
        # sxy = 0.0

        eyy = p_bulk_density * gacc * (locs[:, 1] - ymax) / (p_undrained_lambda + 2 * p_mu)

        syy = p_bulk_density * gacc * (locs[:, 1] - ymax)
        sxx = (p_undrained_lambda + 2.0 * p_mu) * eyy
        sxy = 0.0
        szz = p_undrained_lambda / (2 * p_undrained_lambda + 2 * p_mu) * (sxx + syy)

        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = sxx
        stress[0,:, 1] = syy
        stress[0,:, 2] = szz
        stress[0,:, 3] = sxy
        return stress

    def gacc(self, locs):
        """Compute gravitational acceleration at locations.
        """
        (npts, dim) = locs.shape
        gravacc = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        gravacc[0,:, 1] = -gacc
        return gravacc


# End of file
