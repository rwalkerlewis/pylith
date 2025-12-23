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

"""Unit tests for QuasistaticThermoelasticity scales."""

import unittest

from pylith.scales.General import General
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.QuasistaticThermoelasticity import QuasistaticThermoelasticity

from pythia.pyre.units.length import meter, km
from pythia.pyre.units.time import second
from pythia.pyre.units.pressure import pascal
from pythia.pyre.units.mass import kg
from pythia.pyre.units.temperature import kelvin
from pythia.pyre.units.power import watt
from pythia.pyre.units.energy import joule


def _get_value(v):
    """Extract numeric value from Pyre unit object or return float directly."""
    return v.value if hasattr(v, 'value') else v


class TestThermoelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for thermoelasticity."""

    def test_setQuasistaticThermoelasticity(self):
        """Test setQuasistaticThermoelasticity with default parameters."""
        scales = General()
        scales._configure()
        
        # Default parameters (as floats for this test)
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
        
        # Check scales  
        self.assertAlmostEqual(lengthScale, _get_value(scales.getLengthScale()), places=5)
        
        # Verify time scale
        expectedTime = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        self.assertAlmostEqual(expectedTime, _get_value(scales.getTimeScale()), places=5)

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
        
        self.assertAlmostEqual(lengthScale, _get_value(scales.getLengthScale()), places=5)
        
        # Verify thermal diffusion time scale
        expectedTime = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        self.assertAlmostEqual(expectedTime, _get_value(scales.getTimeScale()), places=5)

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
        self.assertAlmostEqual(expected, _get_value(timeScale), places=5)
        
        # For typical crustal values, this should be on the order of 10^15-10^16 seconds (~30 million years)
        self.assertTrue(_get_value(timeScale) > 1.0e+13)
        self.assertTrue(_get_value(timeScale) < 1.0e+17)

    def test_getTemperatureScale(self):
        """Test getTemperatureScale."""
        scales = General()
        scales._configure()
        temperatureScale = 100.0 * kelvin
        scales.setTemperatureScale(temperatureScale)
        
        result = ElasticityScales.getTemperatureScale(scales)
        self.assertAlmostEqual(_get_value(temperatureScale), _get_value(result), places=10)

    def test_getHeatFluxScale(self):
        """Test getHeatFluxScale."""
        scales = General()
        scales._configure()
        
        # Set up typical scales with units
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
        self.assertTrue(_get_value(heatFluxScale) > 0.0)
        
        # For typical crustal values, should be reasonable (order of W/m^2)
        # q = k * T / L where k ~ rho * c * L^2 / t
        density = (_get_value(rigidityScale) * _get_value(timeScale) * _get_value(timeScale)) / (_get_value(lengthScale) * _get_value(lengthScale))
        thermalConductivity = (density * _get_value(lengthScale) * _get_value(lengthScale)) / _get_value(timeScale)
        expectedFlux = (thermalConductivity * _get_value(temperatureScale)) / _get_value(lengthScale)
        
        self.assertAlmostEqual(expectedFlux, _get_value(heatFluxScale), places=5)

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
        self.assertTrue(_get_value(stressScale) > 0.0)
        self.assertTrue(_get_value(strainScale) > 0.0)
        self.assertTrue(_get_value(temperatureScale) > 0.0)
        self.assertTrue(_get_value(heatFluxScale) > 0.0)
        
        # Verify stress scale is reasonable for crustal values
        # stress ~ rigidity * displacement / length
        expectedStress = _get_value(scales.getRigidityScale()) * _get_value(scales.getDisplacementScale()) / _get_value(scales.getLengthScale())
        self.assertAlmostEqual(expectedStress, _get_value(stressScale), places=5)
        
        # Verify strain scale
        expectedStrain = _get_value(scales.getDisplacementScale()) / _get_value(scales.getLengthScale())
        self.assertAlmostEqual(expectedStrain, _get_value(strainScale), places=10)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestThermoelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
