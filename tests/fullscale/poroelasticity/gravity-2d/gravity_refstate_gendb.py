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
# @file tests/fullscale/linearporoelasticity/gravity-2d/gravity_refstate_gendb.py
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
        y = numpy.arange(-4000.0, 4000.1, 8000.0)
        x = numpy.zeros(y.shape)
        npts = y.shape[0]

        xy = numpy.vstack((x, y)).transpose()

        from gravity_refstate_soln import AnalyticalSoln
        from gravity_refstate_soln import p_solid_density, p_fluid_density, p_fluid_viscosity, p_porosity, p_shear_modulus, p_drained_bulk_modulus, p_biot_coefficient, p_fluid_bulk_modulus, p_solid_bulk_modulus, p_isotropic_permeability
        soln = AnalyticalSoln()
        stress = soln.stress(xy)
        strain = soln.strain(xy)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 1,
            'values': [
                {
                    'name': "solid_density",
                    'units': "kg/m**3",
                    'data': p_solid_density * numpy.ones((npts,)),
                }, {
                    'name': "fluid_density",
                    'units': "kg/m**3",
                    'data': p_fluid_density * numpy.ones((npts,)),
                }, {
                    'name': "fluid_viscosity",
                    'units': "Pa*s",
                    'data': p_fluid_viscosity * numpy.ones((npts,)),
                }, {
                    'name': "porosity",
                    'units': "none",
                    'data': p_porosity * numpy.ones((npts,)),
                }, {
                    'name': "shear_modulus",
                    'units': "Pa",
                    'data': p_shear_modulus * numpy.ones((npts,)),
                }, {
                    'name': "drained_bulk_modulus",
                    'units': "Pa",
                    'data': p_drained_bulk_modulus * numpy.ones((npts,)),
                }, {
                    'name': "biot_coefficient",
                    'units': "none",
                    'data': p_biot_coefficient * numpy.ones((npts,)),
                }, {
                    'name': "fluid_bulk_modulus",
                    'units': "Pa",
                    'data': p_fluid_bulk_modulus * numpy.ones((npts,)),
                }, {
                    'name': "solid_bulk_modulus",
                    'units': "Pa",
                    'data': p_solid_bulk_modulus * numpy.ones((npts,)),
                }, {
                    'name': "isotropic_permeability",
                    'units': "m**2",
                    'data': p_isotropic_permeability * numpy.ones((npts,)),
                }, {
                    'name': "reference_stress_xx",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 0]),
                }, {
                    'name': "reference_stress_yy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 1]),
                }, {
                    'name': "reference_stress_zz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 2]),
                }, {
                    'name': "reference_stress_xy",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 3]),
                }, {
                    'name': "reference_strain_xx",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 0]),
                }, {
                    'name': "reference_strain_yy",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 1]),
                }, {
                    'name': "reference_strain_zz",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 2]),
                }, {
                    'name': "reference_strain_xy",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 3]),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("gravity_refstate_matfields.spatialdb")
        io.write(data)
        return


# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
