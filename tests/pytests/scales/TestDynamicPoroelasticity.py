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

"""Unit tests for DynamicPoroelasticity scales."""

import unittest

from pylith.scales.scales import Scales
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.DynamicPoroelasticity import DynamicPoroelasticity


class TestDynamicPoroelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for dynamic poroelasticity."""

    def test_setDynamicPoroelasticity_defaults(self):
        """Test setDynamicPoroelasticity with default parameters."""
        scales = Scales()
        
        lengthScale = 100.0e+3  # 100 km
        velocityScale = 3.0e+3  # 3 km/s
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        
        ElasticityScales.setDynamicPoroelasticity(
            scales, lengthScale, velocityScale, permeability, viscosity, rigidity
        )
        
        self.assertAlmostEqual(lengthScale, scales.getLengthScale(), places=5)
        self.assertAlmostEqual(1.0, scales.getDisplacementScale(), places=10)
        self.assertAlmostEqual(rigidity, scales.getRigidityScale(), places=5)
        
        # Time scale based on wave propagation (not diffusion)
        expectedTime = lengthScale / velocityScale
        self.assertAlmostEqual(expectedTime, scales.getTimeScale(), places=5)

    def test_setDynamicPoroelasticity_custom(self):
        """Test setDynamicPoroelasticity with custom parameters."""
        scales = Scales()
        
        lengthScale = 10.0e+3  # 10 km
        velocityScale = 1.5e+3  # 1.5 km/s (sediments)
        permeability = 1.0e-13  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 5.0e+9  # Pa (soft sediments)
        
        ElasticityScales.setDynamicPoroelasticity(
            scales, lengthScale, velocityScale, permeability, viscosity, rigidity
        )
        
        # Time scale from wave propagation
        expectedTime = lengthScale / velocityScale
        self.assertAlmostEqual(expectedTime, scales.getTimeScale(), places=5)
        
        # Should be ~6-7 seconds
        self.assertTrue(expectedTime > 5.0)
        self.assertTrue(expectedTime < 10.0)

    def test_DynamicPoroelasticity_class(self):
        """Test DynamicPoroelasticity convenience class."""
        normalizer = DynamicPoroelasticity()
        normalizer._configure()
        
        # Check scales were set
        self.assertTrue(normalizer.getLengthScale() > 0.0)
        self.assertTrue(normalizer.getTimeScale() > 0.0)
        self.assertTrue(normalizer.getRigidityScale() > 0.0)
        
        # Get poroelastic scales
        pressureScale = ElasticityScales.getFluidPressureScale(normalizer)
        viscosityScale = ElasticityScales.getViscosityScale(normalizer)
        permeabilityScale = ElasticityScales.getPermeabilityScale(normalizer)
        velocityScale = ElasticityScales.getVelocityScale(normalizer)
        
        self.assertTrue(pressureScale > 0.0)
        self.assertTrue(viscosityScale > 0.0)
        self.assertTrue(permeabilityScale > 0.0)
        self.assertTrue(velocityScale > 0.0)

    def test_time_scale_difference_quasi_vs_dynamic(self):
        """Test that dynamic uses wave time, not diffusion time."""
        scalesQuasi = Scales()
        scalesDynamic = Scales()
        
        lengthScale = 10.0e+3  # 10 km
        velocityScale = 3.0e+3  # 3 km/s
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        
        # Quasi-static poroelasticity (diffusion-controlled)
        ElasticityScales.setQuasistaticPoroelasticity(
            scalesQuasi, lengthScale, permeability, viscosity, rigidity
        )
        
        # Dynamic poroelasticity (wave-controlled)
        ElasticityScales.setDynamicPoroelasticity(
            scalesDynamic, lengthScale, velocityScale, permeability, viscosity, rigidity
        )
        
        timeQuasi = scalesQuasi.getTimeScale()
        timeDynamic = scalesDynamic.getTimeScale()
        
        # Dynamic time should be wave travel time
        expectedDynamic = lengthScale / velocityScale
        self.assertAlmostEqual(expectedDynamic, timeDynamic, places=5)
        
        # Quasi-static time should be diffusion time
        expectedQuasi = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        self.assertAlmostEqual(expectedQuasi, timeQuasi, places=5)
        
        # Diffusion time should be much longer than wave time
        self.assertTrue(timeQuasi > timeDynamic * 1000.0)

    def test_integration_dynamic_poroelasticity_workflow(self):
        """Test complete workflow for dynamic poroelasticity."""
        scales = Scales()
        
        lengthScale = 5.0e+3  # 5 km sedimentary basin
        velocityScale = 1.0e+3  # 1 km/s
        permeability = 1.0e-13  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 10.0e+9  # Pa
        
        ElasticityScales.setDynamicPoroelasticity(
            scales, lengthScale, velocityScale, permeability, viscosity, rigidity
        )
        
        # Get all relevant scales
        stressScale = ElasticityScales.getStressScale(scales)
        pressureScale = ElasticityScales.getFluidPressureScale(scales)
        velocityScaleComputed = ElasticityScales.getVelocityScale(scales)
        accelerationScale = ElasticityScales.getAccelerationScale(scales)
        viscosityScale = ElasticityScales.getViscosityScale(scales)
        permeabilityScale = ElasticityScales.getPermeabilityScale(scales)
        densityScale = ElasticityScales.getDensityScale(scales)
        
        # All should be positive
        self.assertTrue(stressScale > 0.0)
        self.assertTrue(pressureScale > 0.0)
        self.assertTrue(velocityScaleComputed > 0.0)
        self.assertTrue(accelerationScale > 0.0)
        self.assertTrue(viscosityScale > 0.0)
        self.assertTrue(permeabilityScale > 0.0)
        self.assertTrue(densityScale > 0.0)
        
        # Pressure should equal stress
        self.assertAlmostEqual(stressScale, pressureScale, places=10)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDynamicPoroelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
