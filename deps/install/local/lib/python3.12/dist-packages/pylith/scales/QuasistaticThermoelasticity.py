# =================================================================================================
# This code is part of SpatialData, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/spatialdata).
#
# Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .General import General


class QuasistaticThermoelasticity(General):
    """
    Convenience object for nondimensionalizing quasi-static thermoelasticity problems.

    Implements `General`.
    """

    DOC_CONFIG = {
        "cfg": """
            [normalizer]
            length_scale = 100.0*km
            displacement_scale = 1.0*m
            shear_modulus = 25.0*GPa
            thermal_conductivity = 2.5*W/(m*K)
            density = 2500.0*kg/m**3
            specific_heat = 1000.0*J/(kg*K)
            temperature_scale = 1.0*K
            """,
    }

    import pythia.pyre.inventory

    from pythia.pyre.units.pressure import pascal, GPa
    from pythia.pyre.units.length import meter, km
    from pythia.pyre.units.time import year, second
    from pythia.pyre.units.temperature import kelvin

    lengthScale = pythia.pyre.inventory.dimensional("length_scale", default=100.0 * km)
    lengthScale.validator = pythia.pyre.inventory.greater(0.0 * meter)
    lengthScale.meta["tip"] = (
        "Length scale in boundary value problem (size of feature controlling displacement, fault)."
    )

    displacementScale = pythia.pyre.inventory.dimensional(
        "displacement_scale", default=1.0 * meter
    )
    displacementScale.validator = pythia.pyre.inventory.greater(0.0 * meter)
    displacementScale.meta["tip"] = (
        "Nominal displacement scale in boundary value problem."
    )

    shearModulus = pythia.pyre.inventory.dimensional(
        "shear_modulus", default=25.0 * GPa
    )
    shearModulus.validator = pythia.pyre.inventory.greater(0.0 * pascal)
    shearModulus.meta["tip"] = "Nominal shear modulus in boundary value problem."

    thermalConductivity = pythia.pyre.inventory.float(
        "thermal_conductivity", default=2.5
    )
    thermalConductivity.validator = pythia.pyre.inventory.greater(0.0)
    thermalConductivity.meta["tip"] = "Nominal thermal conductivity in W/(m*K) in boundary value problem."

    density = pythia.pyre.inventory.float(
        "density", default=2500.0
    )
    density.validator = pythia.pyre.inventory.greater(0.0)
    density.meta["tip"] = "Nominal density in kg/m^3 in boundary value problem."

    specificHeat = pythia.pyre.inventory.float(
        "specific_heat", default=1000.0
    )
    specificHeat.validator = pythia.pyre.inventory.greater(0.0)
    specificHeat.meta["tip"] = "Nominal specific heat capacity in J/(kg*K) in boundary value problem."

    temperatureScale = pythia.pyre.inventory.dimensional(
        "temperature_scale", default=1.0 * kelvin
    )
    temperatureScale.validator = pythia.pyre.inventory.greater(0.0 * kelvin)
    temperatureScale.meta["tip"] = "Temperature scale for nondimensionalization."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="quasistaticthermoelasticity"):
        """
        Constructor.
        """
        General.__init__(self, name)

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        from .ElasticityScales import ElasticityScales

        General._configure(self)

        ElasticityScales.setQuasistaticThermoelasticity(
            self,
            lengthScale=self.inventory.lengthScale,
            thermalConductivity=self.inventory.thermalConductivity,
            density=self.inventory.density,
            specificHeat=self.inventory.specificHeat,
        )

        self.setDisplacementScale(self.inventory.displacementScale)
        self.setTemperatureScale(self.inventory.temperatureScale)


# FACTORIES ////////////////////////////////////////////////////////////


def normalizer():
    """
    Factory associated with QuasistaticThermoelasticity.
    """
    return QuasistaticThermoelasticity()


# End of file
