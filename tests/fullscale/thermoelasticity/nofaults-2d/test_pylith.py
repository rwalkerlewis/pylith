#!/usr/bin/env python3
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""Test driver for thermoelasticity full-scale tests.
"""

import unittest

from pylith.testing.FullTestApp import FullTestCase, TestAbstract, SuiteVerbosity

import TestThermoBar


def test_suite():
    """Create test suite."""
    suite = unittest.TestSuite()
    
    # Add thermoelastic bar tests
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestThermoBar.TestTri))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestThermoBar.TestQuad))
    
    return suite


class ThermoelasticityApp(TestAbstract):
    """Driver application for running thermoelasticity tests."""

    def __init__(self):
        """Constructor."""
        TestAbstract.__init__(self)
        self.testSuite = test_suite()
        return


# Run test suite if executed as script
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbosity", action="store", dest="verbosity", default=2, type=int)
    args = parser.parse_args()
    
    runner = unittest.TextTestRunner(verbosity=args.verbosity)
    runner.run(test_suite())


# End of file
