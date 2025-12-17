# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .scales import ElasticityScales as ModuleElasticityScales

from pythia.pyre.units.length import km, meter
from pythia.pyre.units.time import year, second
from pythia.pyre.units.pressure import pascal
from pythia.pyre.units.mass import kg
from pythia.pyre.units.temperature import kelvin
from pythia.pyre.units.power import watt
from pythia.pyre.units.energy import joule
from pythia.pyre.units.unit import one


class ElasticityScales(ModuleElasticityScales):
    """
    Nondimensionalization for elasticity related boundary value problems.
    """

    @staticmethod
    def setQuasistaticElasticity(scales, lengthScale=100.0 * km, timeScale=year):
        ModuleElasticityScales.setQuasistaticElasticity(
            scales, lengthScale.value, timeScale.value
        )

    @staticmethod
    def setDynamicElasticity(
        scales, lengthScale=100.0 * km, velocityScale=3.0 * km / second
    ):
        ModuleElasticityScales.setDynamicElasticity(
            scales, lengthScale.value, velocityScale.value
        )

    @staticmethod
    def setDynamicPoroelasticity(
        scales,
        lengthScale=100.0 * km,
        velocityScale=3.0 * km / second,
        permeability=1.0e-12 * meter**2,
        viscosity=1.0e-3 * pascal * second,
        rigidity=25.0e9 * pascal,
    ):
        ModuleElasticityScales.setDynamicPoroelasticity(
            scales,
            lengthScale.value,
            velocityScale.value,
            permeability.value,
            viscosity.value,
            rigidity.value,
        )

    @staticmethod
    def setQuasistaticPoroelasticity(
        scales,
        lengthScale=100.0 * km,
        permeability=1.0e-12 * meter**2,
        viscosity=1.0e-3 * pascal * second,
        rigidity=25.0e9 * pascal,
    ):
        ModuleElasticityScales.setQuasistaticPoroelasticity(
            scales,
            lengthScale.value,
            permeability.value,
            viscosity.value,
            rigidity.value,
        )

    @staticmethod
    def computePoroelasticityTimeScale(viscosity, permeability, length, rigidity):
        timeScale = ModuleElasticityScales.computePoroelasticityTimeScale(
            viscosity.value, permeability.value, length.value, rigidity.value
        )
        return timeScale * second

    @staticmethod
    def setQuasistaticThermoelasticity(
        scales,
        lengthScale=100.0 * km,
        thermalConductivity=2.5 * watt / (meter * kelvin),
        density=2500.0 * kg / meter**3,
        specificHeat=1000.0 * joule / (kg * kelvin),
    ):
        ModuleElasticityScales.setQuasistaticThermoelasticity(
            scales,
            lengthScale.value,
            thermalConductivity.value,
            density.value,
            specificHeat.value,
        )

    @staticmethod
    def computeThermoelasticityTimeScale(
        lengthScale, thermalConductivity, density, specificHeat
    ):
        timeScale = ModuleElasticityScales.computeThermoelasticityTimeScale(
            lengthScale.value,
            thermalConductivity.value,
            density.value,
            specificHeat.value,
        )
        return timeScale * second

    @staticmethod
    def setQuasistaticThermoporoelasticity(
        scales,
        lengthScale=100.0 * km,
        permeability=1.0e-12 * meter**2,
        viscosity=1.0e-3 * pascal * second,
        rigidity=25.0e9 * pascal,
        thermalConductivity=2.5 * watt / (meter * kelvin),
        density=2500.0 * kg / meter**3,
        specificHeat=1000.0 * joule / (kg * kelvin),
    ):
        ModuleElasticityScales.setQuasistaticThermoporoelasticity(
            scales,
            lengthScale.value,
            permeability.value,
            viscosity.value,
            rigidity.value,
            thermalConductivity.value,
            density.value,
            specificHeat.value,
        )

    @staticmethod
    def setHeat(
        scales,
        lengthScale=100.0 * km,
        thermalConductivity=2.5 * watt / (meter * kelvin),
        density=2500.0 * kg / meter**3,
        specificHeat=1000.0 * joule / (kg * kelvin),
    ):
        ModuleElasticityScales.setHeat(
            scales,
            lengthScale.value,
            thermalConductivity.value,
            density.value,
            specificHeat.value,
        )

    @staticmethod
    def computeThermoporoelasticityTimeScale(
        lengthScale,
        permeability,
        viscosity,
        rigidity,
        thermalConductivity,
        density,
        specificHeat,
    ):
        timeScale = ModuleElasticityScales.computeThermoporoelasticityTimeScale(
            lengthScale.value,
            permeability.value,
            viscosity.value,
            rigidity.value,
            thermalConductivity.value,
            density.value,
            specificHeat.value,
        )
        return timeScale * second

    @staticmethod
    def getStressScale(scales):
        return ModuleElasticityScales.getStressScale(scales) * pascal

    @staticmethod
    def getFluidPressureScale(scales):
        return ModuleElasticityScales.getFluidPressureScale(scales) * pascal

    @staticmethod
    def getStrainScale(scales):
        return ModuleElasticityScales.getStrainScale(scales) * one

    @staticmethod
    def getBodyForceScale(scales):
        return ModuleElasticityScales.getBodyForceScale(scales) * pascal / meter

    @staticmethod
    def getDensityScale(scales):
        return ModuleElasticityScales.getDensityScale(scales) * kg / meter**3

    @staticmethod
    def getVelocityScale(scales):
        return ModuleElasticityScales.getVelocityScale(scales) * meter / second

    @staticmethod
    def getAccelerationScale(scales):
        return ModuleElasticityScales.getAccelerationScale(scales) * meter / second**2

    @staticmethod
    def getViscosityScale(scales):
        return ModuleElasticityScales.getViscosityScale(scales) * pascal * second

    @staticmethod
    def getPermeabilityScale(scales):
        return ModuleElasticityScales.getPermeabilityScale(scales) * meter**2

    @staticmethod
    def getTemperatureScale(scales):
        return ModuleElasticityScales.getTemperatureScale(scales) * kelvin

    @staticmethod
    def getHeatFluxScale(scales):
        return ModuleElasticityScales.getHeatFluxScale(scales) * watt / meter**2


# End of file
