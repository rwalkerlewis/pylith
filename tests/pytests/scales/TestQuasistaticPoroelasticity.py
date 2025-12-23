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

"""Unit tests for QuasistaticPoroelasticity scales."""

import unittest

from pylith.scales.scales import Scales
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.QuasistaticPoroelasticity import QuasistaticPoroelasticity


class TestQuasistaticPoroelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for quasi-static poroelasticity."""

    def test_setQuasistaticPoroelasticity_defaults(self):
        """Test setQuasistaticPoroelasticity with default parameters."""
        scales = Scales()
        
        lengthScale = 100.0e+3  # 100 km
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        
        ElasticityScales.setQuasistaticPoroelasticity(
            scales, lengthScale, permeability, viscosity, rigidity
        )
        
        self.assertAlmostEqual(lengthScale, scales.getLengthScale(), places=5)
        self.assertAlmostEqual(1.0, scales.getDisplacementScale(), places=10)
        self.assertAlmostEqual(rigidity, scales.getRigidityScale(), places=5)
        
        # Check time scale (fluid diffusion)
        expectedTime = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        self.assertAlmostEqual(expectedTime, scales.getTimeScale(), places=5)

    def test_computePoroelasticityTimeScale(self):
        """Test computePoroelasticityTimeScale."""
        viscosity = 1.0e-3  # Pa*s
        permeability = 1.0e-12  # m^2
        lengthScale = 100.0e+3  # 100 km
        rigidity = 25.0e+9  # Pa
        
        timeScale = ElasticityScales.computePoroelasticityTimeScale(
            viscosity, permeability, lengthScale, rigidity
        )
        
        expected = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        self.assertAlmostEqual(expected, timeScale, places=5)
        
        # Should be on order of years to thousands of years
        self.assertTrue(timeScale > 1.0e+8)  # > ~3 years
        self.assertTrue(timeScale < 1.0e+15)  # < ~30 million years

    def test_QuasistaticPoroelasticity_class(self):
        """Test QuasistaticPoroelasticity convenience class."""
        normalizer = QuasistaticPoroelasticity()
        normalizer._configure()
        
        # Check scales
        self.assertTrue(normalizer.getLengthScale() > 0.0)
        self.assertTrue(normalizer.getTimeScale() > 0.0)
        
        # Get derived scales
        pressureScale = ElasticityScales.getFluidPressureScale(normalizer)
        viscosityScale = ElasticityScales.getViscosityScale(normalizer)
        permeabilityScale = ElasticityScales.getPermeabilityScale(normalizer)
        
        self.assertTrue(pressureScale > 0.0)
        self.assertTrue(viscosityScale > 0.0)
        self.assertTrue(permeabilityScale > 0.0)

    def test_integration_poroelasticity_workflow(self):
        """Test complete workflow for poroelasticity."""
        scales = Scales()
        
        lengthScale = 50.0e+3  # 50 km
        permeability = 1.0e-13  # m^2 (low permeability)
        viscosity = 1.0e-3  # Pa*s
        rigidity = 30.0e+9  # Pa
        
        ElasticityScales.setQuasistaticPoroelasticity(
            scales, lengthScale, permeability, viscosity, rigidity
        )
        
        # Get various scales
        stressScale = ElasticityScales.getStressScale(scales)
        pressureScale = ElasticityScales.getFluidPressureScale(scales)
        strainScale = ElasticityScales.getStrainScale(scales)
        viscosityScale = ElasticityScales.getViscosityScale(scales)
        permeabilityScale = ElasticityScales.getPermeabilityScale(scales)
        
        # All positive
        self.assertTrue(stressScale > 0.0)
        self.assertTrue(pressureScale > 0.0)
        self.assertTrue(strainScale > 0.0)
        self.assertTrue(viscosityScale > 0.0)
        self.assertTrue(permeabilityScale > 0.0)
        
        # Pressure should equal stress scale
        self.assertAlmostEqual(stressScale, pressureScale, places=10)
        
        # Permeability scale should be L^2
        expectedPerm = lengthScale * lengthScale
        self.assertAlmostEqual(expectedPerm, permeabilityScale, places=5)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestQuasistaticPoroelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
