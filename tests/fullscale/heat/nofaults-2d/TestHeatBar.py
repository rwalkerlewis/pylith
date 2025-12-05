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
# @file tests/fullscale/heat/nofaults-2d/TestHeatBar.py
#
# @brief Test suite for testing pylith with steady-state heat conduction in a bar.

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check

import meshes
import heatbar_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": heatbar_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["temperature"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["heat_material"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=[
                    "density",
                    "specific_heat",
                    "thermal_conductivity",
                ],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestTri(TestCase):

    def setUp(self):
        self.name = "heatbar_tri"
        self.mesh = meshes.Tri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["heatbar_tri.cfg"])


# -------------------------------------------------------------------------------------------------
class TestQuad(TestCase):

    def setUp(self):
        self.name = "heatbar_quad"
        self.mesh = meshes.Quad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["heatbar_quad.cfg"])


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestTri,
        TestQuad,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
