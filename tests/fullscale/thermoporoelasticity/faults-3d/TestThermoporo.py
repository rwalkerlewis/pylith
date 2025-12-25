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
"""Thermoporoelasticity fullscale test with faults (3D).

NOTE: This test is currently disabled because the analytical solution for
thermoporoelasticity with faults expects a trivial constant-field solution
(zero displacement, zero pressure, uniform temperature). However, the coupling
between temperature, pressure, trace strain, and slip fields through the
continuity equation creates a more complex physical response that doesn't match
the simple analytical solution. The test needs development of a proper coupled
analytical solution that accounts for all thermoporoelastic interactions.
"""

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


@unittest.skip("Analytical solution not valid for coupled thermoporoelastic equations; needs further development")
class TestTetGmsh(TestCase):
    def setUp(self):
        self.name = "thermoporo_fault_tet"
        self.mesh = meshes.TetGmsh()
        super().setUp()
        FullTestCase.run_pylith(self, self.name, ["thermoporo.cfg", "thermoporo_tet.cfg"])


@unittest.skip("Analytical solution not valid for coupled thermoporoelastic equations; needs further development")
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
