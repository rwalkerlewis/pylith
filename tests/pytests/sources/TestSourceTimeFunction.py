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

from pylith.testing.TestCases import TestAbstractComponent, TestComponent, make_suite, configureComponent
import pylith.sources.SourceTimeFunction as stf


class TestSourceTimeFunction(TestAbstractComponent):
    """Unit testing of SourceTimeFunction abstract base class."""
    _class = stf.SourceTimeFunction


class TestSourceTimeStep(TestComponent):
    """Unit testing of SourceTimeStep object."""
    _class = stf.SourceTimeStep
    _factory = stf.source_time_function_step

    def test_evaluate_before_origin(self):
        """Test evaluate() before origin time returns 0."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        
        result = obj.evaluate(0.5)
        self.assertEqual(0.0, result)

    def test_evaluate_after_origin(self):
        """Test evaluate() after origin time returns magnitude."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj.magnitude = 2.5
        
        result = obj.evaluate(1.5)
        self.assertEqual(2.5, result)

    def test_evaluate_at_origin(self):
        """Test evaluate() at origin time returns magnitude."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj.magnitude = 1.0
        
        result = obj.evaluate(1.0)
        self.assertEqual(1.0, result)


class TestSourceTimeRamp(TestComponent):
    """Unit testing of SourceTimeRamp object."""
    _class = stf.SourceTimeRamp
    _factory = stf.source_time_function_ramp

    def test_evaluate_before_origin(self):
        """Test evaluate() before origin time returns 0."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._riseTimeValue = 2.0
        obj.magnitude = 1.0
        
        result = obj.evaluate(0.5)
        self.assertEqual(0.0, result)

    def test_evaluate_during_ramp(self):
        """Test evaluate() during ramp returns interpolated value."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._riseTimeValue = 2.0
        obj.magnitude = 4.0
        
        # At t=2.0, which is 1.0 seconds into the 2.0 second ramp
        result = obj.evaluate(2.0)
        self.assertAlmostEqual(2.0, result, places=10)

    def test_evaluate_after_ramp(self):
        """Test evaluate() after ramp returns full magnitude."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._riseTimeValue = 2.0
        obj.magnitude = 4.0
        
        result = obj.evaluate(5.0)
        self.assertEqual(4.0, result)

    def test_derivative_during_ramp(self):
        """Test evaluateDerivative() during ramp returns constant rate."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._riseTimeValue = 2.0
        obj.magnitude = 4.0
        
        result = obj.evaluateDerivative(2.0)
        self.assertAlmostEqual(2.0, result, places=10)  # 4.0 / 2.0


class TestSourceTimeGaussian(TestComponent):
    """Unit testing of SourceTimeGaussian object."""
    _class = stf.SourceTimeGaussian
    _factory = stf.source_time_function_gaussian

    def test_evaluate_before_origin(self):
        """Test evaluate() before origin time returns 0."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._centerTimeValue = 0.5
        obj._sigmaValue = 0.1
        obj.magnitude = 1.0
        
        result = obj.evaluate(0.5)
        self.assertEqual(0.0, result)

    def test_evaluate_at_peak(self):
        """Test evaluate() at peak returns magnitude."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 0.0
        obj._centerTimeValue = 1.0
        obj._sigmaValue = 0.1
        obj.magnitude = 2.0
        
        result = obj.evaluate(1.0)
        self.assertAlmostEqual(2.0, result, places=10)

    def test_evaluate_symmetry(self):
        """Test Gaussian is symmetric about center."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 0.0
        obj._centerTimeValue = 1.0
        obj._sigmaValue = 0.2
        obj.magnitude = 1.0
        
        result_before = obj.evaluate(0.8)
        result_after = obj.evaluate(1.2)
        self.assertAlmostEqual(result_before, result_after, places=10)


class TestSourceTimeRicker(TestComponent):
    """Unit testing of SourceTimeRicker object."""
    _class = stf.SourceTimeRicker
    _factory = stf.source_time_function_ricker

    def test_evaluate_before_origin(self):
        """Test evaluate() before origin time returns 0."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 1.0
        obj._peakFrequencyValue = 1.0
        obj._delayValue = 0.5
        obj.magnitude = 1.0
        
        result = obj.evaluate(0.5)
        self.assertEqual(0.0, result)

    def test_evaluate_at_peak(self):
        """Test evaluate() at peak time returns magnitude."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 0.0
        obj._peakFrequencyValue = 1.0
        obj._delayValue = 1.0
        obj.magnitude = 1.0
        
        # Peak of Ricker wavelet is at t = delay
        result = obj.evaluate(1.0)
        self.assertAlmostEqual(1.0, result, places=10)

    def test_evaluate_zero_crossing(self):
        """Test Ricker wavelet has zero crossing at expected location."""
        obj = self._class()
        configureComponent(obj)
        obj._originTimeValue = 0.0
        obj._peakFrequencyValue = 1.0
        obj._delayValue = 1.0
        obj.magnitude = 1.0
        
        # Ricker wavelet crosses zero at t = delay +/- sqrt(1.5) / (pi * f)
        # For f=1, this is approximately +/- 0.39
        t_zero = 1.0 + np.sqrt(1.5) / np.pi
        result = obj.evaluate(t_zero)
        self.assertAlmostEqual(0.0, result, places=5)


def load_tests(loader, tests, pattern):
    TEST_CLASSES = [
        TestSourceTimeFunction,
        TestSourceTimeStep,
        TestSourceTimeRamp,
        TestSourceTimeGaussian,
        TestSourceTimeRicker,
    ]
    return make_suite(TEST_CLASSES, loader)


if __name__ == "__main__":
    unittest.main(verbosity=2)


# End of file
