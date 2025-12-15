# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Analytical solution for thermoporoelasticity gravity problem with temperature gradient.

This problem considers a 2D domain with:
- Gravitational body forces
- Hydrostatic pore pressure
- Linear temperature gradient with depth
- Thermal expansion effects

Domain geometry:
  - Width: 8000 m (x: -4000 to +4000 m)
  - Height: 8000 m (y: -8000 to 0 m)
  - Surface at y = 0 m

Boundary conditions:
  - Ux = 0 on x = ±4000 m (roller boundaries)
  - Uy = 0 on y = -8000 m (fixed base)
  - P = 0 on y = 0 m (drained surface)
  - T = T_ref on y = 0 m (surface temperature)
  - T = T_ref + dT_dy * H on y = -8000 m (bottom temperature)

The analytical solution for steady-state equilibrium:
  - Pressure: p(y) = ρ_f * g * (y_max - y)
  - Temperature: T(y) = T_ref + dT_dy * (y_max - y)
  - Displacement and stress account for both poroelastic and thermal effects
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
y_max = 0.0  # m
y_min = -8000.0  # m
x_min = -4000.0  # m
x_max = 4000.0  # m

# Temperature gradient (temperature increases with depth)
dT_dy = -0.025  # K/m (25°C per km, negative because T increases as y decreases)


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution for thermoporoelasticity gravity problem with temperature gradient."""

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "trace_strain": self.trace_strain,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "porosity": self.porosity,
            "specific_heat": self.specific_heat,
            "thermal_conductivity": self.thermal_conductivity,
            "reference_temperature": self.reference_temperature,
            "thermal_expansion_coefficient": self.thermal_expansion_coefficient,
            "fluid_thermal_expansion": self.fluid_thermal_expansion,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": {
                "bc_disp_xneg": self.displacement_zero,
                "bc_disp_xpos": self.displacement_zero,
                "bc_disp_yneg": self.displacement_zero,
                "bc_press_ypos": self.pressure_zero,
                "bc_temp_ypos": self.temperature_surface,
            },
        }

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations.

        For a thermoporoelastic medium under gravity with thermal effects:
        - The y-displacement integrates the strain from the surface
        - Thermal expansion adds to the mechanical strain
        """
        (npts, _) = locs.shape
        y = locs[:, 1]

        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        # Effective stress coefficient for vertical displacement
        # Under uniaxial strain (ε_xx = 0), the y-displacement comes from:
        # ε_yy = (σ_yy + α*p) / (λ + 2μ) + α_T * (T - T_ref)
        # where σ_yy = -ρ_bulk * g * (y_max - y)

        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

        # Mechanical contribution from gravity and pore pressure
        # u_y = ∫[y to y_max] ε_yy dy
        uy_mech = (
            -0.5
            / (p_lambda + 2 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * ((y_max - y_min) ** 2 - (y_max - y) ** 2)
        )

        # Thermal contribution: ε_thermal = α_T * (T - T_ref)
        # T - T_ref = dT_dy * (y_max - y)
        # u_y_thermal = ∫[y to y_max] α_T * dT_dy * (y_max - y') dy'
        #             = α_T * dT_dy * 0.5 * [(y_max - y_min)² - (y_max - y)²]
        # Note: dT_dy is negative, so thermal expansion causes upward displacement
        uy_thermal = (
            0.5
            * alpha_T
            * (-dT_dy)  # dT_dy is negative, so -dT_dy is positive
            * ((y_max - y_min) ** 2 - (y_max - y) ** 2)
        )

        disp[0, :, 0] = 0.0
        disp[0, :, 1] = uy_mech + uy_thermal
        return disp

    def displacement_zero(self, locs):
        (npts, _) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations (hydrostatic)."""
        (npts, _) = locs.shape
        y = locs[:, 1]

        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        pressure[0, :, 0] = p_fluid_density * g_acc * (y_max - y)
        return pressure

    def pressure_zero(self, locs):
        (npts, _) = locs.shape
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return pressure

    def temperature(self, locs):
        """Compute temperature field at locations (linear gradient)."""
        (npts, _) = locs.shape
        y = locs[:, 1]

        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        # Temperature increases with depth (dT_dy is negative)
        temperature[0, :, 0] = p_reference_temperature + dT_dy * (y_max - y)
        return temperature

    def temperature_surface(self, locs):
        """Compute surface temperature (reference temperature)."""
        (npts, _) = locs.shape
        temperature = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        temperature[0, :, 0] = p_reference_temperature
        return temperature

    def trace_strain(self, locs):
        """Compute trace strain field at locations."""
        (npts, _) = locs.shape
        y = locs[:, 1]
        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        # Volumetric strain = ε_yy (since ε_xx = 0 for uniaxial strain)
        ev_mech = (
            -1.0
            / (p_lambda + 2.0 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * (y_max - y)
        )

        # Thermal contribution to volumetric strain
        ev_thermal = alpha_T * (-dT_dy) * (y_max - y)

        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        trace[0, :, 0] = ev_mech + ev_thermal
        return trace

    def solid_density(self, locs):
        """Compute solid density field at locations."""
        (npts, _) = locs.shape
        density = p_solid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def fluid_density(self, locs):
        """Compute fluid density field at locations."""
        (npts, _) = locs.shape
        density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def fluid_viscosity(self, locs):
        """Compute fluid viscosity field at locations."""
        (npts, _) = locs.shape
        viscosity = p_fluid_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return viscosity

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations."""
        (npts, _) = locs.shape
        shear_modulus = p_shear_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def drained_bulk_modulus(self, locs):
        """Compute drained bulk modulus field at locations."""
        (npts, _) = locs.shape
        bulk_modulus = p_drained_bulk_modulus * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return bulk_modulus

    def biot_coefficient(self, locs):
        """Compute Biot coefficient field at locations."""
        (npts, _) = locs.shape
        biot_coeff = p_biot_coefficient * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coeff

    def biot_modulus(self, locs):
        """Compute Biot modulus field at locations."""
        (npts, _) = locs.shape
        p_solid_bulk_modulus = (
            p_drained_bulk_modulus / (1.0 - p_biot_coefficient)
            if p_biot_coefficient < 1.0
            else 1.0e50
        )
        p_biot_modulus = p_fluid_bulk_modulus / (
            p_porosity
            + (p_biot_coefficient - p_porosity)
            * p_fluid_bulk_modulus
            / p_solid_bulk_modulus
        )
        modulus = p_biot_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return modulus

    def isotropic_permeability(self, locs):
        """Compute permeability field at locations."""
        (npts, _) = locs.shape
        permeability = p_isotropic_permeability * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return permeability

    def porosity(self, locs):
        """Compute porosity field at locations."""
        (npts, _) = locs.shape
        value = p_porosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return value

    def specific_heat(self, locs):
        """Compute specific heat field at locations."""
        (npts, _) = locs.shape
        value = p_specific_heat * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return value

    def thermal_conductivity(self, locs):
        """Compute thermal conductivity field at locations."""
        (npts, _) = locs.shape
        value = p_thermal_conductivity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return value

    def reference_temperature(self, locs):
        """Compute reference temperature field at locations."""
        (npts, _) = locs.shape
        value = p_reference_temperature * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return value

    def thermal_expansion_coefficient(self, locs):
        """Compute thermal expansion coefficient field at locations."""
        (npts, _) = locs.shape
        value = p_thermal_expansion_coefficient * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return value

    def fluid_thermal_expansion(self, locs):
        """Compute fluid thermal expansion coefficient field at locations."""
        (npts, _) = locs.shape
        value = p_fluid_thermal_expansion * numpy.ones(
            (1, npts, 1), dtype=numpy.float64
        )
        return value

    def strain(self, locs):
        """Compute strain field at locations."""
        (npts, _) = locs.shape
        y = locs[:, 1]
        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        eyy_mech = (
            -1
            / (p_lambda + 2 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * (y_max - y)
        )
        eyy_thermal = alpha_T * (-dT_dy) * (y_max - y)
        eyy = eyy_mech + eyy_thermal

        strain[0, :, 0] = 0.0  # ε_xx
        strain[0, :, 1] = eyy  # ε_yy
        strain[0, :, 2] = 0.0  # ε_zz
        strain[0, :, 3] = 0.0  # ε_xy
        return strain

    def stress(self, locs):
        """Compute stress field at locations."""
        (npts, _) = locs.shape
        y = locs[:, 1]
        p_alpha = p_biot_coefficient
        alpha_T = p_thermal_expansion_coefficient

        # Effective stress: σ_ij = σ'_ij - α * p * δ_ij
        # For thermal effects, we also have thermal stress contributions

        syy = -p_bulk_density * g_acc * (y_max - y)

        # For uniaxial strain (ε_xx = 0):
        # σ_xx = λ * ε_yy - α * p + thermal_stress_contribution
        eyy_mech = (
            -1
            / (p_lambda + 2 * p_mu)
            * (p_bulk_density - p_alpha * p_fluid_density)
            * g_acc
            * (y_max - y)
        )
        eyy_thermal = alpha_T * (-dT_dy) * (y_max - y)
        eyy = eyy_mech + eyy_thermal

        p_val = p_fluid_density * g_acc * (y_max - y)
        dT = (-dT_dy) * (y_max - y)  # T - T_ref

        # Thermal stress: For plane strain, thermoelastic stress contribution is:
        # σ_thermal = -3 * K * α_T * (T - T_ref) for hydrostatic case
        # But for uniaxial strain, we need to be careful
        K = p_drained_bulk_modulus
        thermal_stress = -(3 * K * alpha_T * dT) * p_lambda / (3 * K)

        sxx = (
            p_lambda / (p_lambda + 2 * p_mu) * (-p_bulk_density * g_acc * (y_max - y))
            - 2 * p_mu / (p_lambda + 2 * p_mu) * p_alpha * p_fluid_density * g_acc * (y_max - y)
            - (p_lambda + 2 * p_mu / 3) * 3 * alpha_T * dT  # thermal contribution
        )
        szz = sxx  # Plane strain

        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = 0.0
        return stress


# End of file
