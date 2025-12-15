# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Analytical solution for thermoelasticity with prescribed slip and temperature gradient.

This problem considers:
- Two blocks separated by a vertical fault at x = 0
- Prescribed left-lateral slip on the fault
- Linear temperature gradient in x-direction
- Thermal expansion effects

Domain geometry:
  - Width: 8000 m (x: -4000 to +4000 m)
  - Height: 8000 m (y: -4000 to +4000 m)
  - Fault at x = 0 m (vertical)

Boundary conditions:
  - Roller boundaries on x = ±4000 m (ux = 0)
  - Fixed y on y = -4000 m (uy = 0)
  - Free on y = +4000 m
  - T = T_left on x = -4000 m
  - T = T_right on x = +4000 m

Material properties:
  - Density: ρ = 2500 kg/m³
  - Vs = 3000 m/s, Vp = 5196 m/s (Poisson's ratio ≈ 0.25)
  - Thermal conductivity: k = 3.0 W/(m·K)
  - Specific heat: c = 1000 J/(kg·K)
  - Thermal expansion coefficient: α = 1e-5 /K
  - Reference temperature: T_ref = 300 K

The prescribed slip causes a rigid block motion superposed on thermal expansion.
"""

import numpy


# Material properties
DENSITY = 2500.0  # kg/m³
VS = 3000.0  # m/s
VP = 5196.0  # m/s (gives Poisson's ratio ≈ 0.25)
SPECIFIC_HEAT = 1000.0  # J/(kg·K)
THERMAL_CONDUCTIVITY = 3.0  # W/(m·K)
THERMAL_EXPANSION_COEFF = 1.0e-5  # /K
REFERENCE_TEMPERATURE = 300.0  # K

# Elastic moduli
SHEAR_MODULUS = DENSITY * VS * VS  # Pa
BULK_MODULUS = DENSITY * (VP * VP - 4.0/3.0 * VS * VS)  # Pa

# Geometry
x_min = -4000.0  # m
x_max = 4000.0  # m
y_min = -4000.0  # m
y_max = 4000.0  # m
domain_x = x_max - x_min

# Fault slip (left-lateral)
SLIP = 2.0  # m - left-lateral slip

# Temperature boundary conditions
T_LEFT = 350.0  # K
T_RIGHT = 300.0  # K


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution for thermoelasticity with slip and temperature gradient."""

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "temperature": self.temperature,
            "slip": self.slip,
        }

    def getField(self, name, mesh_entity, pts):
        if name in self.fields:
            return self.fields[name](pts)
        raise KeyError(f"Unknown field '{name}'.")

    def displacement(self, locs):
        """Compute displacement field at locations.

        The displacement has two contributions:
        1. Rigid block motion from fault slip
        2. Thermal expansion from temperature gradient

        For left-lateral slip on a vertical fault at x=0:
        - Left block (x < 0): moves in +y direction by slip/2
        - Right block (x > 0): moves in -y direction by slip/2

        For temperature gradient with fixed x-boundaries:
        - Thermal expansion is constrained, resulting in thermal stress
        - But with roller boundaries (only x constrained), y can expand
        """
        (npts, _) = locs.shape
        x = locs[:, 0]

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

        # Fault slip contribution (rigid block motion)
        # Left-lateral slip: +x side moves in -y, -x side moves in +y
        disp[0, :, 1] = numpy.where(x >= 0, -SLIP / 2.0, SLIP / 2.0)

        # Thermal expansion contribution in x-direction
        # With roller boundaries (ux = 0 at x = ±4000), thermal expansion is constrained
        # So ux from thermal expansion is zero
        disp[0, :, 0] = 0.0

        return disp

    def temperature(self, locs):
        """Compute temperature field at locations (linear gradient)."""
        (npts, _) = locs.shape
        x = locs[:, 0]

        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        temperature[0, :, 0] = T_LEFT + (T_RIGHT - T_LEFT) * (x - x_min) / domain_x
        return temperature

    def slip(self, locs):
        """Compute slip field on fault."""
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        # Left-lateral slip (negative of right-lateral, which is first component)
        slip[0, :, 0] = SLIP  # left-lateral component
        slip[0, :, 1] = 0.0  # opening component
        return slip


# End of file
