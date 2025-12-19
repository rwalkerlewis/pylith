# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsThermoelasticity(PetscComponent):
    """
    Auxiliary subfields associated with thermoelasticity material.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_thermoelastic]
            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.specific_heat.basis_order = 0
            auxiliary_subfields.thermal_conductivity.basis_order = 0
            auxiliary_subfields.reference_temperature.basis_order = 0
            auxiliary_subfields.thermal_expansion_coefficient.basis_order = 0
            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.bulk_modulus.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    density = pythia.pyre.inventory.facility("density", family="subfield", factory=Subfield)
    density.meta['tip'] = "Density subfield."

    specificHeat = pythia.pyre.inventory.facility("specific_heat", family="subfield", factory=Subfield)
    specificHeat.meta['tip'] = "Specific heat subfield."

    thermalConductivity = pythia.pyre.inventory.facility("thermal_conductivity", family="subfield", factory=Subfield)
    thermalConductivity.meta['tip'] = "Thermal conductivity subfield."

    referenceTemperature = pythia.pyre.inventory.facility("reference_temperature", family="subfield", factory=Subfield)
    referenceTemperature.meta['tip'] = "Reference temperature subfield."

    thermalExpansionCoeff = pythia.pyre.inventory.facility("thermal_expansion_coefficient", family="subfield", factory=Subfield)
    thermalExpansionCoeff.meta['tip'] = "Thermal expansion coefficient subfield."

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pythia.pyre.inventory.facility("bulk_modulus", family="subfield", factory=Subfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    bodyForce = pythia.pyre.inventory.facility("body_force", family="subfield", factory=Subfield)
    bodyForce.meta['tip'] = "Body force subfield."

    heatSource = pythia.pyre.inventory.facility("heat_source", family="subfield", factory=Subfield)
    heatSource.meta['tip'] = "Heat source subfield."

    gravitationalAcceleration = pythia.pyre.inventory.facility("gravitational_acceleration", family="subfield", factory=Subfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    def __init__(self, name="auxsubfieldsthermoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsThermoelasticity.
    """
    return AuxSubfieldsThermoelasticity()


# End of file
