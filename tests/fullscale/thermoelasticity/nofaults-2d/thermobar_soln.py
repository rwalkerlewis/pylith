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
    
    SPACE_DIM = 2

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
        npts = pts.shape[0]
        # Linear temperature profile from T_LEFT to T_RIGHT
        T = T_LEFT + (T_RIGHT - T_LEFT) * (x - XLEFT) / LENGTH
        temperature = np.zeros((1, npts, 1), dtype=np.float64)
        temperature[0, :, 0] = T
        return temperature

    def displacement(self, pts):
        """Compute displacement field for 2D plane strain thermoelasticity.
        
        For a bar with:
        - Fixed left boundary (u = 0 at x = XLEFT)
        - Free right boundary (zero traction)
        - Free top/bottom boundaries (zero traction)
        - Linear temperature gradient in x
        
        In 2D plane strain:
        - The constraint ε_zz = 0 enhances the effective thermal expansion
        - The effective thermal expansion factor is (1 + ν) for isotropic materials
        - With ν ≈ 0.25 (from our Vs/Vp ratio), (1 + ν) ≈ 1.25
        
        Additionally, due to the fixed left boundary and free y-boundaries,
        there is a non-trivial u_y field that varies with position.
        
        Args:
            pts (numpy.ndarray): Coordinates of points [npts, dim].
            
        Returns:
            numpy.ndarray: Displacement at points [npts, dim].
        """
        x = pts[:, 0]
        y = pts[:, 1]
        npts = pts.shape[0]
        dim = pts.shape[1]
        
        disp = np.zeros((npts, dim), dtype=np.float64)
        
        # Compute Poisson's ratio from wave velocities
        # ν = (Vp² - 2Vs²) / (2(Vp² - Vs²))
        vp2 = VP * VP
        vs2 = VS * VS
        nu = (vp2 - 2.0 * vs2) / (2.0 * (vp2 - vs2))
        
        # Plane strain correction factor for thermal expansion
        # In plane strain with ε_zz = 0, the effective in-plane thermal strain is:
        # ε_eff = (1 + ν) * α * ΔT
        plane_strain_factor = 1.0 + nu
        
        # For linear temperature: T(x) = T_LEFT + slope * (x - XLEFT)
        # where slope = (T_RIGHT - T_LEFT) / LENGTH
        slope = (T_RIGHT - T_LEFT) / LENGTH
        
        # Integral of thermal strain from XLEFT to x:
        # ∫[XLEFT to x] α * (T(x') - T_ref) dx'
        delta_x = x - XLEFT
        T_avg_minus_ref = (T_LEFT - REFERENCE_TEMPERATURE) + slope * delta_x / 2.0
        
        # u_x with plane strain correction
        disp[:, 0] = plane_strain_factor * THERMAL_EXPANSION_COEFF * T_avg_minus_ref * delta_x
        
        # u_y: Due to the thermal expansion and fixed left boundary,
        # there is a y-displacement that varies with x and y.
        # For a first approximation, use plane strain thermal expansion in y:
        # u_y ≈ (1 + ν) * α * (T_local - T_ref) * y
        # where T_local is the local temperature at position x
        T_local = T_LEFT + slope * delta_x
        disp[:, 1] = plane_strain_factor * THERMAL_EXPANSION_COEFF * (T_local - REFERENCE_TEMPERATURE) * y
        
        # However, the fixed boundary at x = XLEFT constrains u_y = 0 there.
        # So we need to subtract the y-displacement at x = XLEFT:
        # u_y(x,y) = u_y_local(x,y) - u_y_local(XLEFT,y) but with smooth transition
        # For simplicity, we scale u_y by (x - XLEFT) / LENGTH to enforce u_y = 0 at left
        disp[:, 1] *= delta_x / LENGTH
        
        disp3 = np.zeros((1, npts, dim), dtype=np.float64)
        disp3[0, :, :] = disp
        return disp3


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
    return np.zeros((1, npts, dim), dtype=np.float64)


def bc_xneg_temperature(pts):
    """Get temperature BC on left boundary."""
    npts = pts.shape[0]
    field = np.zeros((1, npts, 1), dtype=np.float64)
    field[0, :, 0] = T_LEFT
    return field


def bc_xpos_temperature(pts):
    """Get temperature BC on right boundary."""
    npts = pts.shape[0]
    field = np.zeros((1, npts, 1), dtype=np.float64)
    field[0, :, 0] = T_RIGHT
    return field


# End of file
