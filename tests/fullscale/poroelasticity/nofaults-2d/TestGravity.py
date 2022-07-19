#!/usr/bin/env nemesis
#
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

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check)

import meshes
import gravity_soln

# We do not include trace_strain in the check of the solution fields, because of the
# poor convergence of the series solution.
SOLUTION_FIELDS = ["displacement", "pressure", "cauchy_strain", "cauchy_stress"]
SOLUTION_TOLERANCE = 0.2

# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with gravitational body forces (no initial stress).
    """
    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": gravity_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "bc_ypos", "points"],
                vertex_fields=["displacement", "pressure"],
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
            Check(
                mesh_entities=["poroelastic_xpos", "poroelastic_xneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=[
                    "biot_coefficient",
                    "biot_modulus",
                    "drained_bulk_modulus",
                    "fluid_density",
                    "fluid_viscosity",
                    "isotropic_permeability",
                    "porosity",
                    "shear_modulus",
                    "solid_density",
                ],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic_xpos", "poroelastic_xneg"],
                vertex_fields = ["displacement", "pressure", "cauchy_strain", "cauchy_stress"],
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos", "bc_yneg"],
                vertex_fields=["displacement", "pressure"],
                defaults=defaults,
                tolerance=SOLUTION_TOLERANCE,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        self.name = "gravity_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        self.name = "gravity_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuad,
        TestTri,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
