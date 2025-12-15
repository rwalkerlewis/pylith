# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

import numpy


class AnalyticalSoln(object):
    """Trivial analytical solution for thermoporoelasticity (3D, constant fields)."""

    SPACE_DIM = 3

    def __init__(self, pressure=0.0, temperature=300.0):
        self._p0 = pressure
        self._T0 = temperature

    def getField(self, name, mesh_entity, pts):
        (npts, dim) = pts.shape

        if name == "displacement":
            return numpy.zeros((1, npts, dim), dtype=numpy.float64)
        if name == "pressure":
            field = numpy.zeros((1, npts, 1), dtype=numpy.float64)
            field[0, :, 0] = self._p0
            return field
        if name == "temperature":
            field = numpy.zeros((1, npts, 1), dtype=numpy.float64)
            field[0, :, 0] = self._T0
            return field
        if name == "slip":
            return numpy.zeros((1, npts, dim), dtype=numpy.float64)

        raise KeyError(f"Unknown field '{name}'.")


# End of file
