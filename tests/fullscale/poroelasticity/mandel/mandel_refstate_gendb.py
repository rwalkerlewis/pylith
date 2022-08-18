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
# @file tests/fullscale/linearelasticity/nofaults-2d/gravity_refstate_gendb.py
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
        xy_data = numpy.ones((len(x1) * len(y1), 2), dtype=numpy.float64)
        xy[:, 0] = x.ravel()
        xy[:, 1] = y.ravel()

        from mandel_refstate_soln import AnalyticalSoln
        from mandel_refstate_soln import p_solid_density, p_fluid_density, p_fluid_viscosity, p_porosity, p_shear_modulus, p_drained_bulk_modulus, p_biot_coefficient, p_fluid_bulk_modulus, p_solid_bulk_modulus, p_isotropic_permeability
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
        data_mat = {
            'points': xy,
            'coordsys': cs,
            'data_dim': 2,
            'values': [
                {
                    'name': "solid_density",
                    'units': "kg/m**3",
                    'data': numpy.ravel(p_solid_density * ones_scalar),
                }, {
                    'name': "fluid_density",
                    'units': "kg/m**3",
                    'data': numpy.ravel(p_fluid_density * ones_scalar),
                }, {
                    'name': "fluid_viscosity",
                    'units': "Pa*s",
                    'data': numpy.ravel(p_fluid_viscosity * ones_scalar),
                }, {
                    'name': "porosity",
                    'units': "none",
                    'data': numpy.ravel(p_porosity * ones_scalar),
                }, {
                    'name': "shear_modulus",
                    'units': "Pa",
                    'data': numpy.ravel(p_shear_modulus * ones_scalar),
                }, {
                    'name': "drained_bulk_modulus",
                    'units': "Pa",
                    'data': numpy.ravel(p_drained_bulk_modulus * ones_scalar),
                }, {
                    'name': "biot_coefficient",
                    'units': "none",
                    'data': numpy.ravel(p_biot_coefficient * ones_scalar),
                }, {
                    'name': "fluid_bulk_modulus",
                    'units': "Pa",
                    'data': numpy.ravel(p_fluid_bulk_modulus * ones_scalar),
                }, {
                    'name': "solid_bulk_modulus",
                    'units': "Pa",
                    'data': numpy.ravel(p_solid_bulk_modulus * ones_scalar),
                }, {
                    'name': "isotropic_permeability",
                    'units': "m**2",
                    'data': numpy.ravel(p_isotropic_permeability * ones_scalar),
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
        io_mat = createWriter("mandel_refstate_matfields.spatialdb")
        io_mat.write(data_mat)

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
                    'data': numpy.ravel(soln.sigma_zz(xy)),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io_mat = createWriter("mandel_refstate_ypos_neu.spatialdb")
        io_mat.write(data_ypos)

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
        io_ic = createWriter("mandel_refstate_ic.spatialdb")
        io_ic.write(data_ic) 

        return



# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
