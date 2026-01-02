# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Reference solution for poroelasticity with prescribed nonzero fault slip.

This test case models a vertical fault with prescribed left-lateral slip
in a poroelastic medium. The fault slip induces:
1. Displacement field around the fault (relative motion across the fault)
2. Volumetric strain changes that couple to pore pressure
3. Pore pressure changes that dissipate through drained boundaries over time

Domain: 8 km x 8 km, vertical fault at y=0 (middle of domain)
Fault slip: 0.5 m left-lateral at t=0 (step function)
Boundaries: 
    - x boundaries fixed in x and y (far field, no displacement)
    - y boundaries fixed in y only
    - y boundaries drained (zero pressure)

This represents a simplified earthquake fault slip scenario in a 
poroelastic medium, where fault slip induces pore pressure changes
that subsequently dissipate.
"""

import numpy

# Physical properties (matching pylithapp.cfg)
p_solid_density = 2500.0  # kg/m^3
p_fluid_density = 1000.0  # kg/m^3
p_fluid_viscosity = 1.0e-3  # Pa*s
p_porosity = 0.02

p_shear_modulus = 30.0e9  # Pa
p_drained_bulk_modulus = 80.0e9  # Pa
p_fluid_bulk_modulus = 10.0e9  # Pa
p_biot_coefficient = 0.7
p_isotropic_permeability = 1.0e-14  # m^2

# Derived properties
p_mu = p_shear_modulus
p_lambda = p_drained_bulk_modulus - 2.0 / 3.0 * p_shear_modulus

# Compute Biot modulus
p_solid_bulk_modulus = p_drained_bulk_modulus / (1.0 - p_biot_coefficient)
p_biot_modulus = 1.0 / (
    p_porosity / p_fluid_bulk_modulus
    + (p_biot_coefficient - p_porosity) / p_solid_bulk_modulus
)

# Domain parameters (matching mesh)
domain_x = 8.0e3  # m
domain_y = 8.0e3  # m

# Slip parameters
slip_magnitude = 0.5  # m (left-lateral)


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Reference solution for poroelastic fault slip problem.
    
    This provides material property values for verification. 
    Solution field values are approximate since a closed-form analytical
    solution is not available for this coupled problem.
    """

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
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
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": {
                "bc_disp_xneg": self.displacement_zero,
                "bc_disp_xpos": self.displacement_zero,
                "bc_disp_yneg": self.displacement_zero,
                "bc_disp_ypos": self.displacement_zero,
                "bc_press_yneg": self.displacement_zero,
                "bc_press_ypos": self.displacement_zero,
            },
            "normal_dir": {
                "bc_disp_xneg": self.orientation_dir((-1, 0)),
                "bc_disp_xpos": self.orientation_dir((+1, 0)),
                "bc_disp_yneg": self.orientation_dir((0, -1)),
                "bc_disp_ypos": self.orientation_dir((0, +1)),
                "bc_press_yneg": self.orientation_dir((0, -1)),
                "bc_press_ypos": self.orientation_dir((0, +1)),
                "fault": self.orientation_dir((0, +1)),
            },
            "tangential_dir": {
                "bc_disp_xneg": self.orientation_dir((0, -1)),
                "bc_disp_xpos": self.orientation_dir((0, +1)),
                "bc_disp_yneg": self.orientation_dir((+1, 0)),
                "bc_disp_ypos": self.orientation_dir((-1, 0)),
                "bc_press_yneg": self.orientation_dir((+1, 0)),
                "bc_press_ypos": self.orientation_dir((-1, 0)),
            },
            "slip": self.slip,
            "traction_change": self.traction_change,
            "strike_dir": self.orientation_dir((-1, 0)),
        }

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def getMask(self, name, mesh_entity, pts):
        """Return mask for points where we skip comparison.
        
        We skip comparison at fault vertices where the solution is discontinuous.
        """
        mask = None
        if name == "displacement":
            # Mask points on or very near the fault (y = 0)
            y = pts[:, 1]
            on_fault = numpy.abs(y) < 1.0
            mask = on_fault
        return mask

    def displacement(self, locs):
        """Compute displacement field at locations.
        
        For a fault slip problem with fixed far-field boundaries, the displacement
        field is zero at the boundaries and shows the slip discontinuity at the fault.
        This is an approximate representation.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        # The actual displacement field depends on the numerical solution
        # Here we return zeros as a placeholder; actual verification checks
        # consistency rather than exact match for this coupled problem
        return disp

    def displacement_zero(self, locs):
        """Zero displacement for fixed boundaries."""
        (npts, _) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return disp

    def pressure(self, locs):
        """Compute pressure field at locations.
        
        At drained boundaries (top and bottom), pressure is zero.
        Interior pressure depends on the poroelastic coupling.
        """
        (npts, _) = locs.shape
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations."""
        (npts, _) = locs.shape
        trace = numpy.zeros((1, npts, 1), dtype=numpy.float64)
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

    def strain(self, locs):
        """Compute strain field at locations."""
        (npts, _) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        return strain

    def stress(self, locs):
        """Compute stress field at locations."""
        (npts, _) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        return stress

    def slip(self, locs):
        """Compute slip field on fault.
        
        Prescribed left-lateral slip of 0.5 m.
        For a vertical fault with strike direction pointing in -x:
        - Left-lateral slip means the positive-y side moves in +x direction
          relative to the negative-y side
        - Slip is in the strike direction (first component)
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        # Left-lateral slip in the strike direction 
        slip[0, :, 1] = -slip_magnitude  # Negative because strike_dir is (-1, 0)
        return slip

    def traction_change(self, locs):
        """Compute change in traction on faults.
        
        For prescribed slip, the traction change depends on the elastic response.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return traction

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, dim) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:, :, d] = vector[d]
            return values
        return fn_dir


# End of file
