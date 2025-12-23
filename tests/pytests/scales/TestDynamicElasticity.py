#!/usr/bin/env nemesis
#
# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

"""Unit tests for DynamicElasticity scales."""

import unittest

from pylith.scales.General import General
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.DynamicElasticity import DynamicElasticity


class TestDynamicElasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for dynamic elasticity."""

    def test_setDynamicElasticity_defaults(self):
        """Test setDynamicElasticity with default parameters."""
        scales = General()
        scales._configure()
        
        lengthScale = 100.0e+3  # 100 km
        velocityScale = 3.0e+3  # 3 km/s
        
        ElasticityScales.setDynamicElasticity(scales, lengthScale, velocityScale)
        
        # Check scales
        self.assertAlmostEqual(lengthScale, scales.getLengthScale().value, places=5)
        self.assertAlmostEqual(1.0, scales.getDisplacementScale().value, places=10)
        self.assertAlmostEqual(2.25e+10, scales.getRigidityScale().value, places=5)
        
        # Time scale for wave propagation
        expectedTime = lengthScale / velocityScale
        self.assertAlmostEqual(expectedTime, scales.getTimeScale().value, places=5)

    def test_setDynamicElasticity_custom(self):
        """Test setDynamicElasticity with custom parameters."""
        scales = General()
        scales._configure()
        
        lengthScale = 50.0e+3  # 50 km
        velocityScale = 5.0e+3  # 5 km/s (fast velocity)
        
        ElasticityScales.setDynamicElasticity(scales, lengthScale, velocityScale)
        
        self.assertAlmostEqual(lengthScale, scales.getLengthScale().value, places=5)
        
        # Verify wave travel time
        expectedTime = lengthScale / velocityScale
        self.assertAlmostEqual(expectedTime, scales.getTimeScale().value, places=5)
        
        # Time should be ~10 seconds for these parameters
        self.assertTrue(expectedTime > 5.0)
        self.assertTrue(expectedTime < 15.0)

    def test_DynamicElasticity_class(self):
        """Test DynamicElasticity convenience class."""
        normalizer = DynamicElasticity()
        normalizer._configure()
        
        # Check that scales were set
        lengthScale = normalizer.getLengthScale()
        timeScale = normalizer.getTimeScale()
        
        self.assertTrue(lengthScale.value > 0.0)
        self.assertTrue(timeScale.value > 0.0)
        
        # Get derived scales
        stressScale = ElasticityScales.getStressScale(normalizer)
        velocityScale = ElasticityScales.getVelocityScale(normalizer)
        accelerationScale = ElasticityScales.getAccelerationScale(normalizer)
        
        self.assertTrue(stressScale.value > 0.0)
        self.assertTrue(velocityScale.value > 0.0)
        self.assertTrue(accelerationScale.value > 0.0)

    def test_integration_dynamic_workflow(self):
        """Test complete workflow for dynamic elasticity."""
        scales = General()
        scales._configure()
        
        lengthScale = 10.0e+3  # 10 km
        velocityScale = 3.0e+3  # 3 km/s
        
        ElasticityScales.setDynamicElasticity(scales, lengthScale, velocityScale)
        
        # Get various scales
        stressScale = ElasticityScales.getStressScale(scales)
        strainScale = ElasticityScales.getStrainScale(scales)
        velocityScaleComputed = ElasticityScales.getVelocityScale(scales)
        accelerationScale = ElasticityScales.getAccelerationScale(scales)
        densityScale = ElasticityScales.getDensityScale(scales)
        
        # All scales should be positive
        self.assertTrue(stressScale.value > 0.0)
        self.assertTrue(strainScale.value > 0.0)
        self.assertTrue(velocityScaleComputed.value > 0.0)
        self.assertTrue(accelerationScale.value > 0.0)
        self.assertTrue(densityScale.value > 0.0)
        
        # Verify relationships
        timeScale = scales.getTimeScale()
        displacement = scales.getDisplacementScale()
        
        expectedVelocity = displacement / timeScale
        self.assertAlmostEqual(expectedVelocity.value, velocityScaleComputed.value, places=5)
        
        expectedAcceleration = displacement / (timeScale * timeScale)
        self.assertAlmostEqual(expectedAcceleration.value, accelerationScale.value, places=5)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDynamicElasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
