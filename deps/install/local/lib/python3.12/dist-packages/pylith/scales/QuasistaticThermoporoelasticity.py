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


class QuasistaticThermoporoelasticity(General):
    """
    Convenience object for nondimensionalizing quasi-static thermoporoelasticity problems.

    Implements `General`.
    """

    DOC_CONFIG = {
        "cfg": """
            [normalizer]
            length_scale = 100.0*km
            displacement_scale = 50.0*km
            shear_modulus = 25.0*GPa
            viscosity = 0.001*Pa*s
            permeability = 1.0e-12*m**2
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
        "shear_modulus", default=10.0 * GPa
    )
    shearModulus.validator = pythia.pyre.inventory.greater(0.0 * pascal)
    shearModulus.meta["tip"] = "Nominal shear modulus in boundary value problem."

    viscosity = pythia.pyre.inventory.dimensional(
        "viscosity", default=0.001 * pascal * second
    )
    viscosity.validator = pythia.pyre.inventory.greater(0.0 * pascal * second)
    viscosity.meta["tip"] = "Nominal fluid viscosity in boundary value problem."

    permeability = pythia.pyre.inventory.dimensional(
        "permeability", default=1.0e-13 * meter**2
    )
    permeability.validator = pythia.pyre.inventory.greater(0.0 * meter**2)
    permeability.meta["tip"] = "Nominal permeability in boundary value problem."

    # Thermal properties
    from pythia.pyre.units.power import watt

    thermalConductivity = pythia.pyre.inventory.dimensional(
        "thermal_conductivity", default=2.5 * watt / (meter * kelvin)
    )
    thermalConductivity.validator = pythia.pyre.inventory.greater(
        0.0 * watt / (meter * kelvin)
    )
    thermalConductivity.meta["tip"] = (
        "Nominal thermal conductivity in boundary value problem."
    )

    from pythia.pyre.units.mass import kilogram

    density = pythia.pyre.inventory.dimensional(
        "density", default=2500.0 * kilogram / meter**3
    )
    density.validator = pythia.pyre.inventory.greater(0.0 * kilogram / meter**3)
    density.meta["tip"] = "Nominal density in boundary value problem."

    from pythia.pyre.units.energy import joule

    specificHeat = pythia.pyre.inventory.dimensional(
        "specific_heat", default=1000.0 * joule / (kilogram * kelvin)
    )
    specificHeat.validator = pythia.pyre.inventory.greater(
        0.0 * joule / (kilogram * kelvin)
    )
    specificHeat.meta["tip"] = (
        "Nominal specific heat capacity in boundary value problem."
    )

    temperatureScale = pythia.pyre.inventory.dimensional(
        "temperature_scale", default=1.0 * kelvin
    )
    temperatureScale.validator = pythia.pyre.inventory.greater(0.0 * kelvin)
    temperatureScale.meta["tip"] = (
        "Temperature scale in boundary value problem."
    )

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="quasistaticthermoporoelasticity"):
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

        ElasticityScales.setQuasistaticThermoporoelasticity(
            self,
            lengthScale=self.inventory.lengthScale,
            permeability=self.inventory.permeability,
            viscosity=self.inventory.viscosity,
            rigidity=self.inventory.shearModulus,
            thermalConductivity=self.inventory.thermalConductivity,
            density=self.inventory.density,
            specificHeat=self.inventory.specificHeat,
        )

        self.setDisplacementScale(self.inventory.displacementScale)
        self.setTemperatureScale(self.inventory.temperatureScale)


# FACTORIES ////////////////////////////////////////////////////////////


def normalizer():
    """
    Factory associated with QuasistaticThermoporoelasticity.
    """
    return QuasistaticThermoporoelasticity()


# End of file
