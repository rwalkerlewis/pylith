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

"""Unit tests for ElasticityScales thermoelasticity functions."""

import unittest

from pylith.scales.General import General
from pylith.scales.ElasticityScales import ElasticityScales
from pythia.pyre.units.length import meter
from pythia.pyre.units.pressure import pascal
from pythia.pyre.units.time import second
from pythia.pyre.units.temperature import kelvin


class TestThermoelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for thermoelasticity."""

    def test_setQuasistaticThermoelasticity(self):
        """Test setQuasistaticThermoelasticity with default parameters."""
        scales = General()
        scales._configure()
        
        # Default parameters
        lengthScale = 100.0e+3  # 100 km in meters
        thermalConductivity = 2.5  # W/(m*K)
        density = 2500.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        ElasticityScales.setQuasistaticThermoelasticity(
            scales, 
            lengthScale, 
            thermalConductivity, 
            density, 
            specificHeat
        )
        
        # Check that scales were set
        self.assertAlmostEqual(lengthScale, scales.getLengthScale().value, places=5)
        self.assertAlmostEqual(1.0, scales.getDisplacementScale().value, places=10)
        self.assertAlmostEqual(2.5e+10, scales.getRigidityScale().value, places=5)
        self.assertAlmostEqual(1.0, scales.getTemperatureScale().value, places=10)
        
        # Check time scale was computed (thermal diffusion time)
        expectedTime = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        self.assertAlmostEqual(expectedTime, scales.getTimeScale().value, places=5)

    def test_setQuasistaticThermoelasticity_custom(self):
        """Test setQuasistaticThermoelasticity with custom parameters."""
        scales = General()
        scales._configure()
        
        # Custom parameters
        lengthScale = 50.0e+3  # 50 km
        thermalConductivity = 3.0  # W/(m*K)
        density = 3000.0  # kg/m^3
        specificHeat = 800.0  # J/(kg*K)
        
        ElasticityScales.setQuasistaticThermoelasticity(
            scales,
            lengthScale,
            thermalConductivity,
            density,
            specificHeat
        )
        
        self.assertAlmostEqual(lengthScale, scales.getLengthScale().value, places=5)
        
        # Verify thermal diffusion time scale
        expectedTime = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        self.assertAlmostEqual(expectedTime, scales.getTimeScale().value, places=5)

    def test_computeThermoelasticityTimeScale(self):
        """Test computeThermoelasticityTimeScale."""
        lengthScale = 100.0e+3  # 100 km
        thermalConductivity = 2.5  # W/(m*K)
        density = 2500.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        timeScale = ElasticityScales.computeThermoelasticityTimeScale(
            lengthScale,
            thermalConductivity,
            density,
            specificHeat
        )
        
        # Expected: t = rho * c * L^2 / k
        expected = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        self.assertAlmostEqual(expected, timeScale.value, places=5)
        
        # For typical crustal values, this should be on the order of 10^14-10^16 seconds (~millions of years)
        self.assertTrue(timeScale.value > 1.0e+13)
        self.assertTrue(timeScale.value < 1.0e+17)

    def test_getTemperatureScale(self):
        """Test getTemperatureScale."""
        scales = General()
        scales._configure()
        temperatureScale = 100.0 * kelvin
        scales.setTemperatureScale(temperatureScale)
        
        result = ElasticityScales.getTemperatureScale(scales)
        self.assertAlmostEqual(temperatureScale.value, result.value, places=10)

    def test_getHeatFluxScale(self):
        """Test getHeatFluxScale."""
        scales = General()
        scales._configure()
        
        # Set up typical scales
        lengthScale = 100.0e+3 * meter
        temperatureScale = 100.0 * kelvin
        rigidityScale = 2.5e+10 * pascal
        timeScale = 1.0e+14 * second
        displacementScale = 1.0 * meter
        
        scales.setLengthScale(lengthScale)
        scales.setTemperatureScale(temperatureScale)
        scales.setRigidityScale(rigidityScale)
        scales.setTimeScale(timeScale)
        scales.setDisplacementScale(displacementScale)
        
        heatFluxScale = ElasticityScales.getHeatFluxScale(scales)
        
        # Heat flux scale should be positive
        self.assertTrue(heatFluxScale.value > 0.0)
        
        # For typical crustal values, should be reasonable (order of W/m^2)
        # q = k * T / L where k ~ rho * c * L^2 / t
        density = (rigidityScale.value * timeScale.value * timeScale.value) / (lengthScale.value * lengthScale.value)
        thermalConductivity = (density * lengthScale.value * lengthScale.value) / timeScale.value
        expectedFlux = (thermalConductivity * temperatureScale.value) / lengthScale.value
        
        self.assertAlmostEqual(expectedFlux, heatFluxScale.value, places=5)

    def test_integration_thermoelasticity_workflow(self):
        """Test complete workflow for setting up thermoelasticity scales."""
        scales = General()
        scales._configure()
        
        # Set up for a crustal thermoelastic problem
        lengthScale = 100.0e+3  # 100 km
        thermalConductivity = 2.5  # W/(m*K)
        density = 2700.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        # Initialize scales
        ElasticityScales.setQuasistaticThermoelasticity(
            scales,
            lengthScale,
            thermalConductivity,
            density,
            specificHeat
        )
        
        # Get various scales for nondimensionalization
        stressScale = ElasticityScales.getStressScale(scales)
        strainScale = ElasticityScales.getStrainScale(scales)
        temperatureScale = ElasticityScales.getTemperatureScale(scales)
        heatFluxScale = ElasticityScales.getHeatFluxScale(scales)
        
        # Verify all scales are positive
        self.assertTrue(stressScale.value > 0.0)
        self.assertTrue(strainScale.value > 0.0)
        self.assertTrue(temperatureScale.value > 0.0)
        self.assertTrue(heatFluxScale.value > 0.0)
        
        # Verify stress scale is reasonable for crustal values
        # stress ~ rigidity * displacement / length
        expectedStress = scales.getRigidityScale() * scales.getDisplacementScale() / scales.getLengthScale()
        self.assertAlmostEqual(expectedStress.value, stressScale.value, places=5)
        
        # Verify strain scale
        expectedStrain = scales.getDisplacementScale() / scales.getLengthScale()
        self.assertAlmostEqual(expectedStrain, strainScale.value, places=10)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestThermoelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
