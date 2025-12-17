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
# @file tests/pytests/scales/TestHeat.py
#
# @brief Unit testing of Python Heat object.

import unittest


class TestHeat(unittest.TestCase):
    """Unit testing of Heat object.
    """

    def test_constructor(self):
        """Test constructor.
        """
        from pylith.scales.Heat import Heat
        scale = Heat()
        self.assertTrue(scale is not None)

    def test_preinitialize(self):
        """Test preinitialize.
        """
        from pylith.scales.Heat import Heat
        from pylith.testing.UnitTestApp import configureComponent
        from pylith.problems.Problem import Problem

        scale = Heat()
        configureComponent(scale)

        problem = Problem()
        configureComponent(problem)

        scale.preinitialize(problem)

        # Check that scales were set
        normalizer = problem.normalizer
        self.assertAlmostEqual(normalizer.getLengthScale(), 100.0e+3, places=1)
        self.assertAlmostEqual(normalizer.getTemperatureScale(), 1.0, places=10)
        self.assertTrue(normalizer.getTimeScale() > 0.0)

    def test_inventory(self):
        """Test inventory.
        """
        from pylith.scales.Heat import Heat
        from pylith.testing.UnitTestApp import configureComponent
        from pythia.pyre.units.length import km, meter
        from pythia.pyre.units.mass import kg
        from pythia.pyre.units.temperature import kelvin
        from pythia.pyre.units.power import watt
        from pythia.pyre.units.energy import joule

        scale = Heat()
        scale.inventory.lengthScale = 10.0*km
        scale.inventory.thermalConductivity = 3.0*watt/(meter*kelvin)
        scale.inventory.density = 2700.0*kg/meter**3
        scale.inventory.specificHeat = 900.0*joule/(kg*kelvin)
        scale.inventory.temperatureScale = 10.0*kelvin
        configureComponent(scale)

        # Verify inventory was set correctly
        self.assertAlmostEqual(scale.lengthScale.value, 10.0e+3, places=1)
        self.assertAlmostEqual(scale.thermalConductivity.value, 3.0, places=10)
        self.assertAlmostEqual(scale.density.value, 2700.0, places=1)
        self.assertAlmostEqual(scale.specificHeat.value, 900.0, places=1)
        self.assertAlmostEqual(scale.temperatureScale.value, 10.0, places=10)

    def test_time_scale(self):
        """Test time scale calculation.
        """
        from pylith.scales.Heat import Heat
        from pylith.testing.UnitTestApp import configureComponent
        from pylith.problems.Problem import Problem
        from pythia.pyre.units.length import meter
        from pythia.pyre.units.mass import kg
        from pythia.pyre.units.temperature import kelvin
        from pythia.pyre.units.power import watt
        from pythia.pyre.units.energy import joule

        # Set up specific values to compute expected time scale
        # Time scale = rho * c * L^2 / k
        L = 1000.0  # meters
        k = 2.5  # W/(m*K)
        rho = 2500.0  # kg/m^3
        c = 1000.0  # J/(kg*K)
        expected_time = (rho * c * L * L) / k  # seconds

        scale = Heat()
        scale.inventory.lengthScale = L*meter
        scale.inventory.thermalConductivity = k*watt/(meter*kelvin)
        scale.inventory.density = rho*kg/meter**3
        scale.inventory.specificHeat = c*joule/(kg*kelvin)
        configureComponent(scale)

        problem = Problem()
        configureComponent(problem)

        scale.preinitialize(problem)

        # Check that time scale matches expected thermal diffusion time scale
        normalizer = problem.normalizer
        actual_time = normalizer.getTimeScale()
        self.assertAlmostEqual(actual_time, expected_time, places=1)


def test_suite():
    """Test suite.
    """
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestHeat))
    return suite


def test_classes():
    """Test classes.
    """
    return [TestHeat]


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
