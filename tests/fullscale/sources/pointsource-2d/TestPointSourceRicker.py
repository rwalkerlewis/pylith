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

from pylith.testing.FullTestApp import run_pylith
from pylith.testing import has_h5py


class TestPointSourceRicker(unittest.TestCase):
    """Full-scale test for a 2D moment-tensor point source with a Ricker wavelet."""

    def setUp(self):
        run_pylith("pointsource_ricker", ["pylithapp.cfg", "pointsource_ricker.cfg"], nprocs=1)

    def test_nonzero_displacement(self):
        if not has_h5py():
            return

        import numpy
        import h5py

        filename = "output/pointsource_ricker-domain.h5"
        with h5py.File(filename, "r") as h5:
            self.assertTrue("vertex_fields" in h5.keys())
            self.assertTrue("displacement" in h5["vertex_fields"].keys())
            disp = h5["vertex_fields/displacement"][:]
            self.assertGreater(numpy.max(numpy.abs(disp)), 0.0)


def test_cases():
    return [TestPointSourceRicker]


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file

