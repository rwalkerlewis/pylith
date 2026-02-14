#!/usr/bin/env nemesis
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
Analytical solution for point source test.

For an isotropic (explosion) point source in an elastic medium, the displacement
field radiates outward from the source location. This solution provides an
approximation based on the Green's function for a point source.
"""

import numpy as np


# Physical parameters
DENSITY = 2500.0  # kg/m^3
VS = 1000.0  # m/s (shear wave speed)
VP = 1732.0  # m/s (P-wave speed)

# Source parameters
SOURCE_X = 0.0  # m
SOURCE_Y = 0.0  # m
MAGNITUDE = 1.0e+15  # N*m
FREQUENCY = 2.0  # Hz
DELAY = 1.0  # s
ORIGIN_TIME = 0.0  # s


def ricker_wavelet(t, origin_time, frequency, delay):
    """Compute Ricker wavelet at time t."""
    tau = t - origin_time - delay
    pi2_f2_tau2 = (np.pi * frequency * tau) ** 2
    return (1.0 - 2.0 * pi2_f2_tau2) * np.exp(-pi2_f2_tau2)


class AnalyticalSoln:
    """Analytical solution for point source problem."""

    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        """Initialize analytical solution."""
        self.source_location = np.array([SOURCE_X, SOURCE_Y])
        self.magnitude = MAGNITUDE
        self.frequency = FREQUENCY
        self.delay = DELAY
        self.origin_time = ORIGIN_TIME

    def getField(self, name, mesh_entity, pts):
        """Get field values at points.

        Args:
            name (str): Name of field.
            mesh_entity: Mesh entity info.
            pts: NumPy array of coordinates.

        Returns:
            NumPy array of field values.
        """
        if name == "displacement":
            return self.displacement(pts)
        elif name == "velocity":
            return self.velocity(pts)
        else:
            raise ValueError(f"Unknown field '{name}'.")

    def displacement(self, locs, t=5.0):
        """Compute displacement at locations for time t.

        For an explosion source, displacement radiates outward from the source.

        Args:
            locs: NumPy array of coordinates [npts, ndim].
            t: Time (default 5.0 s).

        Returns:
            NumPy array of displacement [npts, ndim].
        """
        npts = locs.shape[0]
        disp = np.zeros((npts, self.SPACE_DIM), dtype=np.float64)

        for i in range(npts):
            x = locs[i, 0]
            y = locs[i, 1]

            # Distance from source
            dx = x - SOURCE_X
            dy = y - SOURCE_Y
            r = np.sqrt(dx**2 + dy**2)

            if r < 1.0e-10:
                continue

            # Travel time for P-wave
            travel_time = r / VP
            retarded_time = t - travel_time

            if retarded_time < ORIGIN_TIME:
                continue

            # Source time function
            stf = ricker_wavelet(retarded_time, ORIGIN_TIME, FREQUENCY, DELAY)

            # Amplitude (simplified Green's function)
            # For 2D, Green's function has 1/sqrt(r) behavior
            amplitude = MAGNITUDE * stf / (4.0 * np.pi * DENSITY * VP**2 * np.sqrt(r + 1.0))

            # Radial displacement
            disp[i, 0] = amplitude * dx / r
            disp[i, 1] = amplitude * dy / r

        return disp

    def velocity(self, locs, t=5.0, dt=1.0e-4):
        """Compute velocity at locations for time t.

        Uses finite difference approximation.

        Args:
            locs: NumPy array of coordinates [npts, ndim].
            t: Time (default 5.0 s).
            dt: Time step for finite difference.

        Returns:
            NumPy array of velocity [npts, ndim].
        """
        disp_plus = self.displacement(locs, t + dt)
        disp_minus = self.displacement(locs, t - dt)
        return (disp_plus - disp_minus) / (2.0 * dt)


# End of file
