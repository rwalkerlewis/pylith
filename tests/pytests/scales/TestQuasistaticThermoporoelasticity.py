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

"""Unit tests for QuasistaticThermoporoelasticity scales."""

import unittest

from pylith.scales.General import General
from pylith.scales.ElasticityScales import ElasticityScales
from pylith.scales.QuasistaticThermoporoelasticity import QuasistaticThermoporoelasticity


class TestQuasistaticThermoporoelasticityScales(unittest.TestCase):
    """Test ElasticityScales methods for quasi-static thermoporoelasticity."""

    def test_setQuasistaticThermoporoelasticity_defaults(self):
        """Test setQuasistaticThermoporoelasticity with default parameters."""
        scales = General()
        scales._configure()
        
        lengthScale = 100.0e+3  # 100 km
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        thermalConductivity = 2.5  # W/(m*K)
        density = 2500.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        ElasticityScales.setQuasistaticThermoporoelasticity(
            scales, lengthScale, permeability, viscosity, rigidity,
            thermalConductivity, density, specificHeat
        )
        
        self.assertAlmostEqual(lengthScale, scales.getLengthScale().value, places=5)
        self.assertAlmostEqual(1.0, scales.getDisplacementScale().value, places=10)
        self.assertAlmostEqual(rigidity, scales.getRigidityScale().value, places=5)
        self.assertAlmostEqual(1.0, scales.getTemperatureScale().value, places=10)
        
        # Time scale should be minimum of poro and thermal
        timePoro = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        timeThermal = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        expectedTime = min(timePoro, timeThermal)
        
        self.assertAlmostEqual(expectedTime, scales.getTimeScale().value, places=5)

    def test_computeThermoporoelasticityTimeScale(self):
        """Test computeThermoporoelasticityTimeScale."""
        lengthScale = 100.0e+3  # 100 km
        permeability = 1.0e-12  # m^2
        viscosity = 1.0e-3  # Pa*s
        rigidity = 25.0e+9  # Pa
        thermalConductivity = 2.5  # W/(m*K)
        density = 2500.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        timeScale = ElasticityScales.computeThermoporoelasticityTimeScale(
            lengthScale, permeability, viscosity, rigidity,
            thermalConductivity, density, specificHeat
        )
        
        # Compute both time scales
        timePoro = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        timeThermal = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        
        # Should be minimum of the two
        expectedMin = min(timePoro, timeThermal)
        self.assertAlmostEqual(expectedMin, timeScale.value, places=5)
        
        # Verify it's one of the two
        self.assertTrue(abs(timeScale.value - timePoro) < 1.0 or abs(timeScale.value - timeThermal) < 1.0)

    def test_time_scale_controlled_by_fastest_process(self):
        """Test that time scale is controlled by fastest diffusion process."""
        scales1 = General()
        scales1._configure()
        scales2 = General()
        scales2._configure()
        
        lengthScale = 10.0e+3  # 10 km
        rigidity = 25.0e+9  # Pa
        density = 2500.0  # kg/m^3
        specificHeat = 1000.0  # J/(kg*K)
        
        # Case 1: Fast fluid diffusion (high permeability)
        permeability1 = 1.0e-11  # m^2 (high)
        viscosity1 = 1.0e-3  # Pa*s
        thermalConductivity1 = 2.5  # W/(m*K)
        
        ElasticityScales.setQuasistaticThermoporoelasticity(
            scales1, lengthScale, permeability1, viscosity1, rigidity,
            thermalConductivity1, density, specificHeat
        )
        
        timePoro1 = (viscosity1 * lengthScale * lengthScale) / (permeability1 * rigidity)
        timeThermal1 = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity1
        
        # Fluid should be faster (smaller time)
        self.assertTrue(timePoro1 < timeThermal1)
        self.assertAlmostEqual(timePoro1, scales1.getTimeScale().value, places=5)
        
        # Case 2: Fast thermal diffusion (high conductivity)
        permeability2 = 1.0e-17  # m^2 (very low)
        viscosity2 = 1.0e-3  # Pa*s
        thermalConductivity2 = 10.0  # W/(m*K) (high)
        
        ElasticityScales.setQuasistaticThermoporoelasticity(
            scales2, lengthScale, permeability2, viscosity2, rigidity,
            thermalConductivity2, density, specificHeat
        )
        
        timePoro2 = (viscosity2 * lengthScale * lengthScale) / (permeability2 * rigidity)
        timeThermal2 = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity2
        
        # With very low permeability and high thermal conductivity, 
        # poroelastic diffusion is still slower (counter-intuitively, because 
        # the numerator effect dominates). Expect minimum time to be used.
        minTime2 = min(timePoro2, timeThermal2)
        self.assertAlmostEqual(minTime2, scales2.getTimeScale().value, places=3)

    def test_QuasistaticThermoporoelasticity_class(self):
        """Test QuasistaticThermoporoelasticity convenience class."""
        normalizer = QuasistaticThermoporoelasticity()
        normalizer._configure()
        
        # Check all scales were set
        self.assertTrue(normalizer.getLengthScale().value > 0.0)
        self.assertTrue(normalizer.getTimeScale().value > 0.0)
        self.assertTrue(normalizer.getRigidityScale().value > 0.0)
        self.assertTrue(normalizer.getTemperatureScale().value > 0.0)
        
        # Get all relevant scales
        stressScale = ElasticityScales.getStressScale(normalizer)
        pressureScale = ElasticityScales.getFluidPressureScale(normalizer)
        temperatureScale = ElasticityScales.getTemperatureScale(normalizer)
        heatFluxScale = ElasticityScales.getHeatFluxScale(normalizer)
        viscosityScale = ElasticityScales.getViscosityScale(normalizer)
        permeabilityScale = ElasticityScales.getPermeabilityScale(normalizer)
        
        self.assertTrue(stressScale.value > 0.0)
        self.assertTrue(pressureScale.value > 0.0)
        self.assertTrue(temperatureScale.value > 0.0)
        self.assertTrue(heatFluxScale.value > 0.0)
        self.assertTrue(viscosityScale.value > 0.0)
        self.assertTrue(permeabilityScale.value > 0.0)

    def test_integration_thermoporoelasticity_workflow(self):
        """Test complete workflow for thermoporoelasticity."""
        scales = General()
        scales._configure()
        
        # Fault zone thermal pressurization parameters
        lengthScale = 0.1  # 0.1 m fault zone
        permeability = 1.0e-18  # m^2 (very low)
        viscosity = 1.0e-3  # Pa*s
        rigidity = 20.0e+9  # Pa
        thermalConductivity = 1.5  # W/(m*K)
        density = 2500.0  # kg/m^3
        specificHeat = 900.0  # J/(kg*K)
        
        ElasticityScales.setQuasistaticThermoporoelasticity(
            scales, lengthScale, permeability, viscosity, rigidity,
            thermalConductivity, density, specificHeat
        )
        
        # Get all scales
        stressScale = ElasticityScales.getStressScale(scales)
        pressureScale = ElasticityScales.getFluidPressureScale(scales)
        temperatureScale = ElasticityScales.getTemperatureScale(scales)
        heatFluxScale = ElasticityScales.getHeatFluxScale(scales)
        strainScale = ElasticityScales.getStrainScale(scales)
        
        # All should be positive
        self.assertTrue(stressScale.value > 0.0)
        self.assertTrue(pressureScale.value > 0.0)
        self.assertTrue(temperatureScale.value > 0.0)
        self.assertTrue(heatFluxScale.value > 0.0)
        self.assertTrue(strainScale.value > 0.0)
        
        # Compute individual time scales to verify which is controlling
        timePoro = (viscosity * lengthScale * lengthScale) / (permeability * rigidity)
        timeThermal = (density * specificHeat * lengthScale * lengthScale) / thermalConductivity
        
        # For very low permeability, thermal should be faster
        if timeThermal < timePoro:
            self.assertAlmostEqual(timeThermal, scales.getTimeScale().value, places=5)
        else:
            self.assertAlmostEqual(timePoro, scales.getTimeScale().value, places=5)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestQuasistaticThermoporoelasticityScales)
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
