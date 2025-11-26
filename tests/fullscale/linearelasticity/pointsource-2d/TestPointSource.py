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

from pylith.testing.FullTestApp import (FullTestCase, Check)

import meshes
import pointsource_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Test suite for testing PyLith with point moment tensor sources.
    """

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": pointsource_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "velocity"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["elastic"],
                vertex_fields=["displacement", "velocity"],
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestExplosionQuad(TestCase):
    """Test explosion source with quad mesh."""

    def setUp(self):
        self.name = "pointsource_explosion_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp.cfg", "pointsource_explosion.cfg", "mesh_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestExplosionTri(TestCase):
    """Test explosion source with tri mesh."""

    def setUp(self):
        self.name = "pointsource_explosion_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp.cfg", "pointsource_explosion.cfg", "mesh_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestStrikeslipQuad(TestCase):
    """Test strike-slip source with quad mesh."""

    def setUp(self):
        self.name = "pointsource_strikeslip_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp.cfg", "pointsource_strikeslip.cfg", "mesh_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestStrikeslipTri(TestCase):
    """Test strike-slip source with tri mesh."""

    def setUp(self):
        self.name = "pointsource_strikeslip_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp.cfg", "pointsource_strikeslip.cfg", "mesh_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestExplosionQuad,
        TestExplosionTri,
        TestStrikeslipQuad,
        TestStrikeslipTri,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
