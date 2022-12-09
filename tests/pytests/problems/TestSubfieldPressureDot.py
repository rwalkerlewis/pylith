#!/usr/bin/env nemesis
#
# ======================================================================
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
# ======================================================================
#
# @file tests/pytests/problems/TestSubfieldPressureDot.py
#
# @brief Unit testing of Python SubfieldPressureDot object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.problems.SubfieldPressureDot import (SubfieldPressureDot, soln_subfield)


class TestSubfieldPressureDot(TestComponent):
    """Unit testing of SubfieldPressureDot object.
    """
    _class = SubfieldPressureDot
    _factory = soln_subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfieldPressureDot))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
