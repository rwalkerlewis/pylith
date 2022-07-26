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
# @file tests/fullscale/poroelasticity/gravity-3d/gravity_gendb.py
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
        z1 = numpy.arange(-9000.0, 0.1, 9000.0)
        x1 = numpy.arange(-9000.0, 0.1, 9000.0)
        y1 = numpy.arange(-9000.0, 0.1, 9000.0)
        x, y, z = numpy.meshgrid(x1, y1, z1)
        xyz = numpy.zeros((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
        xyz_data = numpy.ones((len(x1) * len(y1) * len(z1), 3), dtype=numpy.float64)
        xyz[:, 0] = x.ravel()
        xyz[:, 1] = y.ravel()
        xyz[:, 2] = z.ravel()

        from gravity_soln import AnalyticalSoln
        from gravity_soln import p_solid_density, p_fluid_density, p_fluid_viscosity, p_porosity, p_shear_modulus, p_drained_bulk_modulus, p_biot_coefficient, p_fluid_bulk_modulus, p_solid_bulk_modulus, p_isotropic_permeability
        soln = AnalyticalSoln()
        stress = soln.stress(xyz)
        strain = soln.strain(xyz)
        zero_scalar = soln.zero_scalar(xyz)
        ones_scalar = soln.ones_scalar(xyz)        
        pressure = soln.pressure(xyz)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 3
        cs._configure()
        data = {
            'points': xyz,
            'coordsys': cs,
            'data_dim': 1,
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
                    'name': "reference_stress_yz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 4]),
                }, {
                    'name': "reference_stress_xz",
                    'units': "Pa",
                    'data': numpy.ravel(stress[0, :, 5]),
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
                }, {
                    'name': "reference_strain_yz",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 4]),
                }, {
                    'name': "reference_strain_xz",
                    'units': "none",
                    'data': numpy.ravel(strain[0, :, 5]),
                }
            ]
        }

        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("gravity_matfields.spatialdb")
        io.write(data)

        data_ic = {
            'points': xyz,
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
                    'name': "displacement_z",
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
        io_ic = createWriter("gravity_ic.spatialdb")
        io_ic.write(data_ic) 

        return


# ======================================================================
if __name__ == "__main__":
    app = GenerateDB().run()


# End of file
