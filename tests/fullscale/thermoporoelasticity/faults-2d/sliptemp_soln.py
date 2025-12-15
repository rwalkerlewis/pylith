# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Analytical solution for thermoporoelasticity with prescribed slip and temperature gradient.

This problem considers:
- Two blocks separated by a vertical fault at x = 0
- Prescribed left-lateral slip on the fault
- Linear temperature gradient in x-direction
- Constant pore pressure

Domain geometry:
  - Width: 8000 m (x: -4000 to +4000 m)
  - Height: 8000 m (y: -4000 to +4000 m)
  - Fault at x = 0 m (vertical)

Boundary conditions:
  - Ux = ±slip/2 on fault (left-lateral slip)
  - Roller boundaries on x = ±4000 m
  - Fixed y on y = -4000 m
  - P = 0 (constant)
  - T = T_left on x = -4000 m
  - T = T_right on x = +4000 m

The prescribed slip causes a rigid block motion superposed on thermal expansion.
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
p_specific_heat = 1000.0  # J/(kg·K)
p_thermal_conductivity = 3.0  # W/(m·K)
p_reference_temperature = 300.0  # K
p_thermal_expansion_coefficient = 1.0e-5  # 1/K
p_fluid_thermal_expansion = 3.0e-4  # 1/K

# Derived properties
p_mu = p_shear_modulus
p_lambda = p_drained_bulk_modulus - 2.0 / 3.0 * p_shear_modulus

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
    """Analytical solution for thermoporoelasticity with slip and temperature gradient."""

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "trace_strain": self.trace_strain,
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
        - Left block (x < 0): moves in +y direction
        - Right block (x > 0): moves in -y direction

        For temperature gradient with free boundaries:
        - Thermal strain causes x-displacement
        """
        (npts, _) = locs.shape
        x = locs[:, 0]

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

        # Fault slip contribution (rigid block motion)
        # Left-lateral slip: +x side moves in -y, -x side moves in +y
        disp[0, :, 1] = numpy.where(x >= 0, -SLIP / 2.0, SLIP / 2.0)

        # Thermal expansion contribution
        # Temperature is linear from T_LEFT to T_RIGHT
        # T(x) = T_LEFT + (T_RIGHT - T_LEFT) * (x - x_min) / domain_x
        # Thermal strain: ε_th = α * (T - T_ref)
        # For free expansion: u_x = ∫ α * (T - T_ref) dx from x_min to x

        alpha_T = p_thermal_expansion_coefficient
        slope = (T_RIGHT - T_LEFT) / domain_x

        # For the left block (x < 0):
        delta_x_left = x - x_min
        T_avg_left = (T_LEFT - p_reference_temperature) + slope * delta_x_left / 2.0
        ux_thermal_left = alpha_T * T_avg_left * delta_x_left

        # For the right block (x >= 0):
        # The displacement jumps at the fault, but thermal expansion continues
        # We integrate from x_min to x, accounting for the fault
        ux_thermal_right = alpha_T * (
            (T_LEFT - p_reference_temperature + slope * (0 - x_min) / 2.0) * (0 - x_min)
            + (T_LEFT + slope * (0 - x_min) - p_reference_temperature + slope * x / 2.0) * x
        )

        disp[0, :, 0] = numpy.where(x >= 0, ux_thermal_right, ux_thermal_left)

        return disp

    def pressure(self, locs):
        """Compute pressure field at locations (zero)."""
        (npts, _) = locs.shape
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return pressure

    def temperature(self, locs):
        """Compute temperature field at locations (linear gradient)."""
        (npts, _) = locs.shape
        x = locs[:, 0]

        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        temperature[0, :, 0] = T_LEFT + (T_RIGHT - T_LEFT) * (x - x_min) / domain_x
        return temperature

    def trace_strain(self, locs):
        """Compute trace strain field at locations."""
        (npts, _) = locs.shape
        x = locs[:, 0]

        alpha_T = p_thermal_expansion_coefficient
        slope = (T_RIGHT - T_LEFT) / domain_x

        # T - T_ref at each location
        T = T_LEFT + slope * (x - x_min)
        dT = T - p_reference_temperature

        # Volumetric thermal strain
        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        trace[0, :, 0] = alpha_T * dT
        return trace

    def slip(self, locs):
        """Compute slip field on fault."""
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        # Left-lateral slip (negative of right-lateral, which is first component)
        slip[0, :, 0] = SLIP  # left-lateral component
        slip[0, :, 1] = 0.0  # opening component
        return slip


# End of file
