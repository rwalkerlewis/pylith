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

from pylith.scales.General import General
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.QuasistaticPoroelasticity import QuasistaticPoroelasticity


def _get_value(v):
    """Extract numeric value from Pyre unit object or return float directly."""
    return v.value if hasattr(v, 'value') else v


class TestQuasistaticPoroelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for quasi-static poroelasticity."""

    def test_setQuasistaticPoroelasticity_defaults(self):
        """Test setQuasistaticPoroelasticity with default parameters."""
        scales = General()
        scales._configure()
        
        lengthScale = 100.0e+3  # 100 km
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        
        ElasticityScales.setQuasistaticPoroelasticity(
            scales, lengthScale, permeability, viscosity, rigidity
        )
        
        self.assertAlmostEqual(lengthScale, _get_value(scales.getLengthScale()), places=5)
        self.assertAlmostEqual(1.0, _get_value(scales.getDisplacementScale()), places=10)
        self.assertAlmostEqual(rigidity, _get_value(scales.getRigidityScale()), places=5)
        
        # Check time scale (fluid diffusion)
        expectedTime = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        self.assertAlmostEqual(expectedTime, _get_value(scales.getTimeScale()), places=5)

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
        self.assertAlmostEqual(expected, _get_value(timeScale), places=5)
        
        # Should be on order of years to thousands of years
        self.assertTrue(_get_value(timeScale) > 1.0e+8)  # > ~3 years
        self.assertTrue(_get_value(timeScale) < 1.0e+15)  # < ~30 million years

    def test_QuasistaticPoroelasticity_class(self):
        """Test QuasistaticPoroelasticity convenience class."""
        normalizer = QuasistaticPoroelasticity()
        normalizer._configure()
        
        # Check scales
        self.assertTrue(_get_value(normalizer.getLengthScale()) > 0.0)
        self.assertTrue(_get_value(normalizer.getTimeScale()) > 0.0)
        
        # Get derived scales
        pressureScale = ElasticityScales.getFluidPressureScale(normalizer)
        viscosityScale = ElasticityScales.getViscosityScale(normalizer)
        permeabilityScale = ElasticityScales.getPermeabilityScale(normalizer)
        
        self.assertTrue(_get_value(pressureScale) > 0.0)
        self.assertTrue(_get_value(viscosityScale) > 0.0)
        self.assertTrue(_get_value(permeabilityScale) > 0.0)

    def test_integration_poroelasticity_workflow(self):
        """Test complete workflow for poroelasticity."""
        scales = General()
        scales._configure()
        
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
        self.assertTrue(_get_value(stressScale) > 0.0)
        self.assertTrue(_get_value(pressureScale) > 0.0)
        self.assertTrue(_get_value(strainScale) > 0.0)
        self.assertTrue(_get_value(viscosityScale) > 0.0)
        self.assertTrue(_get_value(permeabilityScale) > 0.0)
        
        # Pressure should equal stress scale
        self.assertAlmostEqual(_get_value(stressScale), _get_value(pressureScale), places=10)
        
        # Permeability scale should be L^2
        expectedPerm = lengthScale * lengthScale
        self.assertAlmostEqual(expectedPerm, _get_value(permeabilityScale), places=5)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestQuasistaticPoroelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
