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
import numpy as np

from pylith.testing.TestCases import TestComponent, make_suite, configureComponent
import pylith.sources.PointMomentTensor as pmt


class TestPointMomentTensor(TestComponent):
    """Unit testing of PointMomentTensor object."""
    _class = pmt.PointMomentTensor
    _factory = pmt.point_source

    def test_default_location(self):
        """Test default location is at origin."""
        obj = self._class()
        configureComponent(obj)
        self.assertEqual([0.0, 0.0, 0.0], obj.location)

    def test_default_moment_tensor(self):
        """Test default moment tensor is isotropic."""
        obj = self._class()
        configureComponent(obj)
        self.assertEqual([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], obj.momentTensor)

    def test_validate_location_2d(self):
        """Test location validation for 2D."""
        result = pmt.validateLocation([1.0, 2.0])
        self.assertEqual([1.0, 2.0], result)

    def test_validate_location_3d(self):
        """Test location validation for 3D."""
        result = pmt.validateLocation([1.0, 2.0, 3.0])
        self.assertEqual([1.0, 2.0, 3.0], result)

    def test_validate_location_invalid_length(self):
        """Test location validation fails for invalid length."""
        with self.assertRaises(ValueError):
            pmt.validateLocation([1.0])
        with self.assertRaises(ValueError):
            pmt.validateLocation([1.0, 2.0, 3.0, 4.0])

    def test_validate_moment_tensor_2d(self):
        """Test moment tensor validation for 2D (3 components)."""
        result = pmt.validateMomentTensor([1.0, -1.0, 0.0])
        self.assertEqual([1.0, -1.0, 0.0], result)

    def test_validate_moment_tensor_3d(self):
        """Test moment tensor validation for 3D (6 components)."""
        result = pmt.validateMomentTensor([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        self.assertEqual([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], result)

    def test_validate_moment_tensor_invalid(self):
        """Test moment tensor validation fails for invalid length."""
        with self.assertRaises(ValueError):
            pmt.validateMomentTensor([1.0, 2.0])
        with self.assertRaises(ValueError):
            pmt.validateMomentTensor([1.0, 2.0, 3.0, 4.0, 5.0])


class TestPointMomentTensorMomentMatrix(unittest.TestCase):
    """Test moment tensor matrix conversion."""

    def test_moment_tensor_full_2d(self):
        """Test getMomentTensorFull() for 2D."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        # Set up normalized moment tensor
        obj._momentTensorNd = np.array([1.0, -1.0, 0.5])
        
        M = obj.getMomentTensorFull(2)
        
        # Expected symmetric matrix
        expected = np.array([
            [1.0, 0.5],
            [0.5, -1.0]
        ])
        np.testing.assert_array_almost_equal(expected, M)

    def test_moment_tensor_full_3d(self):
        """Test getMomentTensorFull() for 3D."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        # Set up normalized moment tensor [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
        obj._momentTensorNd = np.array([1.0, 2.0, 3.0, 0.5, 0.25, 0.1])
        
        M = obj.getMomentTensorFull(3)
        
        # Expected symmetric matrix
        expected = np.array([
            [1.0, 0.5, 0.25],
            [0.5, 2.0, 0.1],
            [0.25, 0.1, 3.0]
        ])
        np.testing.assert_array_almost_equal(expected, M)

    def test_moment_tensor_isotropic(self):
        """Test isotropic moment tensor gives diagonal matrix."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        # Isotropic tensor
        obj._momentTensorNd = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        
        M = obj.getMomentTensorFull(3)
        
        # Should be identity-like
        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        np.testing.assert_array_almost_equal(expected, M)

    def test_moment_tensor_double_couple(self):
        """Test double-couple moment tensor."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        # Strike-slip double couple
        obj._momentTensorNd = np.array([1.0, -1.0, 0.0, 0.0, 0.0, 0.0])
        
        M = obj.getMomentTensorFull(3)
        
        # Should have trace = 0
        trace = np.trace(M)
        self.assertAlmostEqual(0.0, trace, places=10)


class TestPointMomentTensorEquivalentForce(unittest.TestCase):
    """Test equivalent force computation."""

    def test_equivalent_force_at_source(self):
        """Test equivalent force computation at source location."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        obj._locationNd = np.array([0.0, 0.0])
        obj._momentTensorNd = np.array([1.0, 1.0, 0.0])
        obj._magnitudeNd = 1.0
        
        # Mock the time function
        class MockTimeFunction:
            def evaluate(self, t):
                return 1.0
        obj.timeFunction = MockTimeFunction()
        
        # At the source, the force should be computed
        x = np.array([0.0, 0.0])
        force = obj.computeEquivalentForce(x, 0.0, 2, charLength=0.1)
        
        # Force should be a 2D vector
        self.assertEqual(2, len(force))

    def test_equivalent_force_far_from_source(self):
        """Test equivalent force far from source is small."""
        obj = pmt.PointMomentTensor()
        configureComponent(obj)
        
        obj._locationNd = np.array([0.0, 0.0])
        obj._momentTensorNd = np.array([1.0, 1.0, 0.0])
        obj._magnitudeNd = 1.0
        
        # Mock the time function
        class MockTimeFunction:
            def evaluate(self, t):
                return 1.0
        obj.timeFunction = MockTimeFunction()
        
        # Far from source, force should be small
        x = np.array([10.0, 10.0])
        force = obj.computeEquivalentForce(x, 0.0, 2, charLength=0.1)
        
        # Force magnitude should be very small
        force_mag = np.sqrt(np.sum(force**2))
        self.assertLess(force_mag, 1e-10)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestPointMomentTensor,
        TestPointMomentTensorMomentMatrix,
        TestPointMomentTensorEquivalentForce,
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
