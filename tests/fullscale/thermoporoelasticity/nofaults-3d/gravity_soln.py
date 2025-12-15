# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Analytical solution for 3D thermoporoelasticity gravity problem with temperature gradient.

This problem considers a 3D domain with:
- Gravitational body forces
- Hydrostatic pore pressure
- Linear temperature gradient with depth (z-direction)
- Thermal expansion effects

Boundary conditions:
  - Ux = 0 on x = ±boundary
  - Uy = 0 on y = ±boundary
  - Uz = 0 on z = z_min
  - P = 0 on z = z_max (drained surface)
  - T = T_ref on z = z_max (surface temperature)

The analytical solution for steady-state equilibrium:
  - Pressure: p(z) = ρ_f * g * (z_max - z)
  - Temperature: T(z) = T_ref + dT_dz * (z_max - z)
  - Displacement accounts for both poroelastic and thermal effects
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
p_bulk_density = (1.0 - p_porosity) * p_solid_density + p_porosity * p_fluid_density

# Geometry
g_acc = 9.80665  # m/s²
z_max = 0.0  # m (surface)
z_min = -16000.0  # m (base, 16 km depth)

# Temperature gradient (temperature increases with depth)
dT_dz = -0.025  # K/m (25°C per km, negative because T increases as z decreases)


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution for 3D thermoporoelasticity gravity problem with temperature gradient."""

    SPACE_DIM = 3
    TENSOR_SIZE = 6

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "trace_strain": self.trace_strain,
        }

    def getField(self, name, mesh_entity, pts):
        if name in self.fields:
            return self.fields[name](pts)
        raise KeyError(f"Unknown field '{name}'.")

    def displacement(self, locs):
        """Compute displacement field at locations."""
        (npts, dim) = locs.shape
        z = locs[:, 2]

        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        disp = numpy.zeros((1, npts, dim), dtype=numpy.float64)

        # Mechanical contribution from gravity and pore pressure
        uz_mech = (
            -0.5
            / (p_lambda + 2 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * ((z_max - z_min) ** 2 - (z_max - z) ** 2)
        )

        # Thermal contribution
        uz_thermal = (
            0.5
            * alpha_T
            * (-dT_dz)
            * ((z_max - z_min) ** 2 - (z_max - z) ** 2)
        )

        disp[0, :, 0] = 0.0
        disp[0, :, 1] = 0.0
        disp[0, :, 2] = uz_mech + uz_thermal
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations (hydrostatic)."""
        (npts, _) = locs.shape
        z = locs[:, 2]

        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        pressure[0, :, 0] = p_fluid_density * g_acc * (z_max - z)
        return pressure

    def temperature(self, locs):
        """Compute temperature field at locations (linear gradient)."""
        (npts, _) = locs.shape
        z = locs[:, 2]

        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        temperature[0, :, 0] = p_reference_temperature + dT_dz * (z_max - z)
        return temperature

    def trace_strain(self, locs):
        """Compute trace strain field at locations."""
        (npts, _) = locs.shape
        z = locs[:, 2]
        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        ev_mech = (
            -1.0
            / (p_lambda + 2.0 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * (z_max - z)
        )

        ev_thermal = alpha_T * (-dT_dz) * (z_max - z)

        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        trace[0, :, 0] = ev_mech + ev_thermal
        return trace


# End of file
