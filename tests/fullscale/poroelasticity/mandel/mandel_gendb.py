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
# @file tests/fullscale/linearporoelasticity/mandel/mandel_refstate_gendb.py
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
        x1 = numpy.arange(-0.1, 10.1, 0.1)
        y1 = numpy.arange(-0.1, 1.01, 0.1)
        x, y = numpy.meshgrid(x1, y1)

        xy = numpy.zeros((len(x1) * len(y1), 2), dtype=numpy.float64)
        xy[:, 0] = x.ravel()
        xy[:, 1] = y.ravel()

        from mandel_soln import AnalyticalSoln
        from mandel_soln import p_solid_density, p_fluid_density, p_fluid_viscosity, p_porosity, p_shear_modulus, p_drained_bulk_modulus, p_biot_coefficient, p_fluid_bulk_modulus, p_solid_bulk_modulus, p_isotropic_permeability
        soln = AnalyticalSoln()
        stress = soln.stress(xy)
        strain = soln.strain(xy)
        ones_scalar = soln.ones_scalar(xy)
        zero_scalar = soln.zero_scalar(xy)
        # pressure = soln.zero_scalar(xy)
        pressure = soln.initial_pressure(xy)

        # Aux Fields
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        # Y+ Neumann BC
        data_ypos = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [
                {
                    'name': "initial_amplitude_tangential",
                    'units': "Pa",
                    'data': numpy.ravel(zero_scalar),
                }, {
                    'name': "initial_amplitude_normal",
                    'units': "Pa",
                    'data': numpy.ravel(zero_scalar),
                }
            ]
        }
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io_mat = createWriter("mandel_ypos_neu.spatialdb")
        io_mat.write(data_ypos)

        # Initial conditions
        data_ic = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [
                {
                    'name': "displacement_x",
                    'units': "m",
                    'data': numpy.ravel(zero_scalar),
                }, {
                    'name': "displacement_y",
                    'units': "m",
                    'data': numpy.ravel(zero_scalar),
                }, {
                    'name': "pressure",
                    'units': "Pa",
                    'data': numpy.ravel(pressure),
                }, {
                    'name': "trace_strain",
                    'units': "none",
                    'data': numpy.ravel(zero_scalar),
                }
            ]
        }
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io_ic = createWriter("mandel_ic.spatialdb")
        io_ic.write(data_ic) 

        return



# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
