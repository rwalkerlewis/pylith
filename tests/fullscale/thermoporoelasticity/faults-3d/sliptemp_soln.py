# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Analytical solution for 3D thermoporoelasticity with prescribed slip and temperature gradient.

This problem considers:
- A 3D domain with a vertical fault
- Prescribed left-lateral slip on the fault
- Linear temperature gradient in x-direction
- Constant pore pressure

The prescribed slip causes rigid block motion.
"""

import numpy


# Physical properties
p_solid_density = 2500.0  # kg/m³
p_fluid_density = 1000.0  # kg/m³
p_fluid_viscosity = 1.0e-3  # Pa·s
p_porosity = 0.02

p_shear_modulus = 30.0e9  # Pa
p_drained_bulk_modulus = 80.0e9  # Pa
p_fluid_bulk_modulus = 10.0e9  # Pa
p_biot_coefficient = 0.7
p_isotropic_permeability = 1.0e-14  # m²

# Thermal properties
p_reference_temperature = 300.0  # K
p_thermal_expansion_coefficient = 1.0e-5  # 1/K

# Fault slip (left-lateral)
SLIP = 2.0  # m - left-lateral slip

# Temperature boundary conditions
T_LEFT = 350.0  # K
T_RIGHT = 300.0  # K

# Domain geometry (approximate)
x_min = -4000.0  # m
x_max = 4000.0  # m
domain_x = x_max - x_min


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution for 3D thermoporoelasticity with slip and temperature gradient."""

    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "slip": self.slip,
        }

    def getField(self, name, mesh_entity, pts):
        if name in self.fields:
            return self.fields[name](pts)
        raise KeyError(f"Unknown field '{name}'.")

    def displacement(self, locs):
        """Compute displacement field at locations.

        For left-lateral slip on a vertical fault at x=0 (strike in y-direction):
        - Left block (x < 0): moves in +y direction
        - Right block (x > 0): moves in -y direction
        """
        (npts, dim) = locs.shape
        x = locs[:, 0]

        disp = numpy.zeros((1, npts, dim), dtype=numpy.float64)

        # Fault slip contribution (rigid block motion)
        disp[0, :, 0] = 0.0
        disp[0, :, 1] = numpy.where(x >= 0, -SLIP / 2.0, SLIP / 2.0)
        disp[0, :, 2] = 0.0

        return disp

    def pressure(self, locs):
        """Compute pressure field at locations (zero)."""
        (npts, _) = locs.shape
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return pressure

    def temperature(self, locs):
        """Compute temperature field at locations (linear gradient in x)."""
        (npts, _) = locs.shape
        x = locs[:, 0]

        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        temperature[0, :, 0] = T_LEFT + (T_RIGHT - T_LEFT) * (x - x_min) / domain_x
        return temperature

    def slip(self, locs):
        """Compute slip field on fault."""
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        # Left-lateral slip (along strike, y-direction in fault coordinates)
        slip[0, :, 0] = SLIP  # left-lateral component
        slip[0, :, 1] = 0.0  # reverse component
        slip[0, :, 2] = 0.0  # opening component
        return slip


# End of file
