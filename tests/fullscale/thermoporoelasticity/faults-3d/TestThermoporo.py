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

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check

import meshes
import thermoporo_soln


class TestCase(FullTestCase):
    """Thermoporoelasticity fullscale test (3D fault, constant solution)."""

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
                mesh_entities=["thermoporoelastic"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["slip"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]


class TestTetGmsh(TestCase):
    def setUp(self):
        self.name = "thermoporo_fault_tet"
        self.mesh = meshes.TetGmsh()
        super().setUp()
        FullTestCase.run_pylith(self, self.name, ["thermoporo.cfg", "thermoporo_tet.cfg"])


class TestHexGmsh(TestCase):
    def setUp(self):
        self.name = "thermoporo_fault_hex"
        self.mesh = meshes.HexGmsh()
        super().setUp()
        FullTestCase.run_pylith(self, self.name, ["thermoporo.cfg", "thermoporo_hex.cfg"])


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestTetGmsh))
    suite.addTests(loader.loadTestsFromTestCase(TestHexGmsh))
    return suite


if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
