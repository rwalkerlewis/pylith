#!/usr/bin/env nemesis
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/linearelasticity/gravity-2d/gravity_gendb.py
#
# @brief Python script to generate spatial database with auxiliary
# fields for test with gravitational body forces and initial
# stress/strain but no displacement.

import numpy


class GenerateDB(object):
    """Python object to generate spatial database with auxiliary fields for
    test with gravitational body forces and initial stress/strain but no
    displacement.
    """

    def __init__(self):
        """Constructor.
        """
        return

    def run(self):
        """Generate the database.
        """
        # Domain
        y1 = numpy.arange(-4000.0, 4000.1, 100.0)
        x1 = numpy.arange(-4000.0, 4000.1, 100.0)
        x, y = numpy.meshgrid(x1, y1)

        xy = numpy.zeros((len(x1) * len(y1), 2), dtype=numpy.float64)
        xy[:, 0] = x.ravel()
        xy[:, 1] = y.ravel()

        from gravity_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xy) * 0.0
        pres = soln.pressure(xy)
        trace_strain = soln.pressure(xy) * 0.0

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [
                {
                    'name': "displacement_x",
                    'units': "kg/m**3",
                    'data': numpy.ravel(disp[0, :, 0]),
                }, {
                    'name': "displacement_y",
                    'units': "kg/m**3",
                    'data': numpy.ravel(disp[0, :, 1]),
                }, {
                    'name': "pressure",
                    'units': "Pa*s",
                    'data': numpy.ravel(pres[:]),
                }, {
                    'name': "trace_strain",
                    'units': "none",
                    'data': numpy.ravel(trace_strain[:]),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("gravity_ic.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
