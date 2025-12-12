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
#
# @file tests/pytests/sources/TestSquareWavelet.py
#
# @brief Unit testing of Python SquareWavelet source time function.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.SquareWavelet import (SquareWavelet, momenttensorforce_sourcetimefunction)


class TestSquareWavelet(TestComponent):
    """Unit testing of SquareWavelet source time function."""
    _class = SquareWavelet
    _factory = momenttensorforce_sourcetimefunction


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSquareWavelet))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file

