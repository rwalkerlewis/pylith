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
# @file tests/pytests/faults/TestFaultCohesiveKinPoro.py
#
# @brief Unit testing of Python FaultCohesiveKinPoro object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.faults.FaultCohesiveKinPoro import (FaultCohesiveKinPoro, fault)


class TestFaultCohesiveKinPoro(TestComponent):
    """Unit testing of FaultCohesiveKinPoro object.
    """
    _class = FaultCohesiveKinPoro
    _factory = fault


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFaultCohesiveKinPoro))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
