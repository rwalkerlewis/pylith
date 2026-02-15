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
"""Thermoelasticity fullscale test with prescribed slip and temperature gradient."""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check
from pylith.testing import TestCases

import meshes
import slipthermal_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Thermoelasticity fullscale test with prescribed slip and temperature gradient."""

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": slipthermal_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["thermoelastic"],
                vertex_fields=["displacement", "temperature"],
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
class TestTriGmsh(TestCase):
    """Test case for triangular mesh."""

    def setUp(self):
        self.name = "slipthermal_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["slipthermal.cfg", "slipthermal_tri.cfg"])


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):
    """Test case for quadrilateral mesh."""

    def setUp(self):
        self.name = "slipthermal_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["slipthermal.cfg", "slipthermal_quad.cfg"])


# -------------------------------------------------------------------------------------------------
def load_tests(loader, tests, pattern):
    TEST_CLASSES = (
        TestTriGmsh,
        TestQuadGmsh,
    )
    return TestCases.make_suite(test_classes=TEST_CLASSES, loader=loader)


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
