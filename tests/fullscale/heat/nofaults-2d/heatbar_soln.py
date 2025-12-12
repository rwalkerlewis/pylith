# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
# @file tests/fullscale/heat/nofaults-2d/heatbar_soln.py
#
# @brief Analytical solution to 2D steady-state heat conduction in a bar.
#
# Steady-state heat conduction with Dirichlet BCs:
#   T(-L,y) = T_left
#   T(+L,y) = T_right
#
# The analytical solution is:
#   T(x,y) = T_left + (T_right - T_left) * (x + L) / (2*L)
#
# For this test:
#   L = 4.0e+3 m (half-width)
#   T_left = 300.0 K
#   T_right = 400.0 K

import numpy

# Physical properties
rho = 2500.0  # kg/m**3, density
c_p = 900.0   # J/(kg*K), specific heat capacity
k = 3.0       # W/(m*K), thermal conductivity

# Domain
L = 4.0e+3  # m, half-width
T_left = 300.0   # K
T_right = 400.0  # K


class AnalyticalSoln(object):
    """Analytical solution to steady-state heat conduction in a bar."""

    SPACE_DIM = 2

    def __init__(self):
        self.fields = {
            "temperature": self.temperature,
            "density": self.density,
            "specific_heat": self.specific_heat,
            "thermal_conductivity": self.thermal_conductivity,
            "initial_amplitude": {
                "bc_xneg": self.bc_xneg_temperature,
                "bc_xpos": self.bc_xpos_temperature,
            },
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def temperature(self, locs):
        """Compute temperature field at locations."""
        (npts, dim) = locs.shape
        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        # Linear temperature profile: T = T_left + (T_right - T_left) * (x + L) / (2*L)
        temperature[0, :, 0] = T_left + (T_right - T_left) * (x + L) / (2.0 * L)
        return temperature

    def density(self, locs):
        """Compute density field at locations."""
        (npts, dim) = locs.shape
        density = rho * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def specific_heat(self, locs):
        """Compute specific heat field at locations."""
        (npts, dim) = locs.shape
        specific_heat = c_p * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return specific_heat

    def thermal_conductivity(self, locs):
        """Compute thermal conductivity field at locations."""
        (npts, dim) = locs.shape
        thermal_conductivity = k * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return thermal_conductivity

    def bc_xneg_temperature(self, locs):
        """Compute Dirichlet BC on -x boundary."""
        (npts, dim) = locs.shape
        bc = T_left * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bc

    def bc_xpos_temperature(self, locs):
        """Compute Dirichlet BC on +x boundary."""
        (npts, dim) = locs.shape
        bc = T_right * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bc


# End of file
