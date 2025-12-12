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
# @file tests/pytests/sources/TestGaussianWavelet.py
#
# @brief Unit testing of Python GaussianWavelet source time function.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.GaussianWavelet import (GaussianWavelet, momenttensorforce_sourcetimefunction)


class TestGaussianWavelet(TestComponent):
    """Unit testing of GaussianWavelet source time function."""
    _class = GaussianWavelet
    _factory = momenttensorforce_sourcetimefunction


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGaussianWavelet))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file

