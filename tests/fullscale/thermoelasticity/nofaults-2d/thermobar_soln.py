# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""
Analytical solution for thermoelastic bar problem.

A bar is fixed on the left boundary and has zero traction on the right boundary.
Temperature is imposed as a linear gradient from left to right.
The thermal strain causes displacement to the right.

Bar geometry:
  - Length L = 8000 m
  - Height H = 4000 m
  - Left boundary: x = -4000 m
  - Right boundary: x = +4000 m

Boundary conditions:
  - Left (x = -L/2): u_x = 0, u_y = 0, T = T_left
  - Right (x = +L/2): Free surface (zero traction), T = T_right
  - Top/Bottom: Zero traction, no heat flux (adiabatic)

Material properties:
  - Density: ρ = 2500 kg/m³
  - Vs = 3000 m/s, Vp = 5196 m/s (Poisson's ratio ≈ 0.25)
  - Thermal conductivity: k = 3.0 W/(m·K)
  - Specific heat: c = 1000 J/(kg·K)
  - Thermal expansion coefficient: α = 1e-5 /K
  - Reference temperature: T_ref = 300 K

Analytical solution for steady-state:
  - Temperature is linear: T(x) = T_left + (T_right - T_left) * (x + L/2) / L
  - Thermal strain: ε_th = α * (T - T_ref) in all directions
  - For free boundary on right, stress is zero throughout (for 1D problem)
  - Total strain equals thermal strain for stress-free state
"""

import numpy as np


# Constants
XLEFT = -4000.0  # m
XRIGHT = +4000.0  # m
LENGTH = XRIGHT - XLEFT
HEIGHT = 4000.0  # m

# Material properties
DENSITY = 2500.0  # kg/m³
VS = 3000.0  # m/s
VP = 5196.0  # m/s (gives Poisson's ratio ≈ 0.25)
SPECIFIC_HEAT = 1000.0  # J/(kg·K)
THERMAL_CONDUCTIVITY = 3.0  # W/(m·K)
THERMAL_EXPANSION_COEFF = 1.0e-5  # /K
REFERENCE_TEMPERATURE = 300.0  # K

# Boundary temperatures
T_LEFT = 350.0  # K (warmer on left)
T_RIGHT = 300.0  # K (at reference on right)

# Elastic moduli
SHEAR_MODULUS = DENSITY * VS * VS  # Pa
BULK_MODULUS = DENSITY * (VP * VP - 4.0/3.0 * VS * VS)  # Pa

# Lame parameters
LAMBDA = BULK_MODULUS - 2.0/3.0 * SHEAR_MODULUS
MU = SHEAR_MODULUS


class AnalyticalSolution:
    """Analytical solution for thermoelastic bar."""

    def __init__(self):
        """Initialize."""
        self.XLEFT = XLEFT
        self.XRIGHT = XRIGHT
        self.fields = ["displacement", "temperature"]
        return

    def getField(self, name, mesh_entity, pts):
        """Get field values at points.
        
        Args:
            name (str): Name of field.
            mesh_entity (str): Type of mesh entity (vertex, cell).
            pts (numpy.ndarray): Coordinates of points.
            
        Returns:
            numpy.ndarray: Field values at points.
        """
        if name == "displacement":
            return self.displacement(pts)
        elif name == "temperature":
            return self.temperature(pts)
        else:
            raise ValueError(f"Unknown field '{name}'.")

    def temperature(self, pts):
        """Compute temperature field (linear gradient).
        
        Args:
            pts (numpy.ndarray): Coordinates of points [npts, dim].
            
        Returns:
            numpy.ndarray: Temperature at points [npts, 1].
        """
        x = pts[:, 0]
        # Linear temperature profile from T_LEFT to T_RIGHT
        T = T_LEFT + (T_RIGHT - T_LEFT) * (x - XLEFT) / LENGTH
        return T.reshape(-1, 1)

    def displacement(self, pts):
        """Compute displacement field.
        
        For a bar with:
        - Fixed left boundary (u = 0 at x = XLEFT)
        - Free right boundary (zero traction)
        - Linear temperature gradient
        
        The analytical solution is derived from:
        - Stress = 0 (free expansion)
        - ε_total = ε_thermal = α * (T - T_ref)
        - u_x(x) = ∫[XLEFT to x] α * (T(x') - T_ref) dx'
        
        Args:
            pts (numpy.ndarray): Coordinates of points [npts, dim].
            
        Returns:
            numpy.ndarray: Displacement at points [npts, dim].
        """
        x = pts[:, 0]
        npts = pts.shape[0]
        dim = pts.shape[1]
        
        disp = np.zeros((npts, dim))
        
        # For linear temperature: T(x) = T_LEFT + slope * (x - XLEFT)
        # where slope = (T_RIGHT - T_LEFT) / LENGTH
        slope = (T_RIGHT - T_LEFT) / LENGTH
        
        # Integral of thermal strain from XLEFT to x:
        # ∫[XLEFT to x] α * (T(x') - T_ref) dx'
        # = α * ∫[XLEFT to x] (T_LEFT - T_ref + slope * (x' - XLEFT)) dx'
        # = α * [(T_LEFT - T_ref) * (x - XLEFT) + slope/2 * (x - XLEFT)²]
        
        delta_x = x - XLEFT
        T_avg_minus_ref = (T_LEFT - REFERENCE_TEMPERATURE) + slope * delta_x / 2.0
        
        disp[:, 0] = THERMAL_EXPANSION_COEFF * T_avg_minus_ref * delta_x
        disp[:, 1] = 0.0  # No y-displacement for 1D thermal expansion with free boundaries
        
        return disp


# Create solution instance for use by test framework
soln = AnalyticalSolution()


def displacement(pts):
    """Get displacement field."""
    return soln.displacement(pts)


def temperature(pts):
    """Get temperature field."""
    return soln.temperature(pts)


# Auxiliary field functions for database
def density(pts):
    """Get density."""
    return DENSITY * np.ones((pts.shape[0], 1))


def vs(pts):
    """Get shear wave velocity."""
    return VS * np.ones((pts.shape[0], 1))


def vp(pts):
    """Get compressional wave velocity."""
    return VP * np.ones((pts.shape[0], 1))


def specific_heat(pts):
    """Get specific heat."""
    return SPECIFIC_HEAT * np.ones((pts.shape[0], 1))


def thermal_conductivity(pts):
    """Get thermal conductivity."""
    return THERMAL_CONDUCTIVITY * np.ones((pts.shape[0], 1))


def reference_temperature(pts):
    """Get reference temperature."""
    return REFERENCE_TEMPERATURE * np.ones((pts.shape[0], 1))


def thermal_expansion_coefficient(pts):
    """Get thermal expansion coefficient."""
    return THERMAL_EXPANSION_COEFF * np.ones((pts.shape[0], 1))


def bc_xneg_displacement(pts):
    """Get displacement BC on left boundary."""
    npts = pts.shape[0]
    dim = pts.shape[1]
    return np.zeros((npts, dim))


def bc_xneg_temperature(pts):
    """Get temperature BC on left boundary."""
    return T_LEFT * np.ones((pts.shape[0], 1))


def bc_xpos_temperature(pts):
    """Get temperature BC on right boundary."""
    return T_RIGHT * np.ones((pts.shape[0], 1))


# End of file
