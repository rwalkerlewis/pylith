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
"""Thermoporoelasticity fullscale test with prescribed slip and temperature gradient (3D)."""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check
from pylith.testing import TestCases

import meshes
import sliptemp_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Thermoporoelasticity fullscale test with prescribed slip and temperature gradient (3D)."""

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": sliptemp_soln.AnalyticalSoln(),
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

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestTetGmsh(TestCase):
    """Test case for tetrahedral mesh."""

    def setUp(self):
        self.name = "sliptemp_tet"
        self.mesh = meshes.TetGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["sliptemp.cfg", "sliptemp_tet.cfg"])


# -------------------------------------------------------------------------------------------------
class TestHexGmsh(TestCase):
    """Test case for hexahedral mesh."""

    def setUp(self):
        self.name = "sliptemp_hex"
        self.mesh = meshes.HexGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["sliptemp.cfg", "sliptemp_hex.cfg"])


# -------------------------------------------------------------------------------------------------
def load_tests(loader, tests, pattern):
    TEST_CLASSES = (
        TestTetGmsh,
        TestHexGmsh,
    )
    return TestCases.make_suite(test_classes=TEST_CLASSES, loader=loader)


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
