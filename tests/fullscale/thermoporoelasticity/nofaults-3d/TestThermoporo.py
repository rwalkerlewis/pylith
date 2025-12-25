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
"""Thermoporoelasticity fullscale test (3D, constant solution).

NOTE: This test is currently disabled due to configuration file parsing issues
with the multi-line array syntax and variable substitution in the 3D test setup.
The configuration needs to be simplified to work with the config parser.
"""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check

import meshes
import thermoporo_soln


class TestCase(FullTestCase):
    """Thermoporoelasticity fullscale test (3D, constant solution)."""

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": thermoporo_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["upper_crust"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["lower_crust"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]


@unittest.skip("Configuration file parsing issues with 3D multi-layer setup; needs further development")
class TestTet(TestCase):
    def setUp(self):
        self.name = "thermoporo_tet"
        self.mesh = meshes.Tet()
        super().setUp()
        FullTestCase.run_pylith(self, self.name, ["thermoporo.cfg", "thermoporo_tet.cfg"])


@unittest.skip("Configuration file parsing issues with 3D multi-layer setup; needs further development")
class TestHex(TestCase):
    def setUp(self):
        self.name = "thermoporo_hex"
        self.mesh = meshes.Hex()
        super().setUp()
        FullTestCase.run_pylith(self, self.name, ["thermoporo.cfg", "thermoporo_hex.cfg"])


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestTet))
    suite.addTests(loader.loadTestsFromTestCase(TestHex))
    return suite


if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
