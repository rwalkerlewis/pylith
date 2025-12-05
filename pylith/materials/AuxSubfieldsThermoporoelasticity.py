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


class AuxSubfieldsThermoporoelasticity(PetscComponent):
    """
    Auxiliary subfields for isotropic, linear thermoporoelasticity.
    
    Combines poroelastic and thermal properties.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat.bulk_rheology.auxiliary_subfields]
            biot_coefficient.basis_order = 0
            biot_modulus.basis_order = 0
            drained_bulk_modulus.basis_order = 0
            shear_modulus.basis_order = 0
            isotropic_permeability.basis_order = 0
            reference_temperature.basis_order = 0
            thermal_expansion_coefficient.basis_order = 0
            fluid_thermal_expansion.basis_order = 0
            thermal_conductivity.basis_order = 0
            specific_heat.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    # Poroelastic properties
    biotCoefficient = pythia.pyre.inventory.facility("biot_coefficient", family="auxiliary_subfield", factory=Subfield)
    biotCoefficient.meta['tip'] = "Biot coefficient subfield."

    biotModulus = pythia.pyre.inventory.facility("biot_modulus", family="auxiliary_subfield", factory=Subfield)
    biotModulus.meta['tip'] = "Biot modulus subfield."

    drainedBulkModulus = pythia.pyre.inventory.facility("drained_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    drainedBulkModulus.meta['tip'] = "Drained bulk modulus subfield."

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    isotropicPermeability = pythia.pyre.inventory.facility("isotropic_permeability", family="auxiliary_subfield", factory=Subfield)
    isotropicPermeability.meta['tip'] = "Isotropic permeability subfield."

    # Thermal properties
    referenceTemperature = pythia.pyre.inventory.facility("reference_temperature", family="auxiliary_subfield", factory=Subfield)
    referenceTemperature.meta['tip'] = "Reference temperature subfield."

    thermalExpansionCoefficient = pythia.pyre.inventory.facility("thermal_expansion_coefficient", family="auxiliary_subfield", factory=Subfield)
    thermalExpansionCoefficient.meta['tip'] = "Thermal expansion coefficient (solid) subfield."

    fluidThermalExpansion = pythia.pyre.inventory.facility("fluid_thermal_expansion", family="auxiliary_subfield", factory=Subfield)
    fluidThermalExpansion.meta['tip'] = "Fluid thermal expansion coefficient subfield."

    thermalConductivity = pythia.pyre.inventory.facility("thermal_conductivity", family="auxiliary_subfield", factory=Subfield)
    thermalConductivity.meta['tip'] = "Thermal conductivity subfield."

    specificHeat = pythia.pyre.inventory.facility("specific_heat", family="auxiliary_subfield", factory=Subfield)
    specificHeat.meta['tip'] = "Specific heat subfield."

    # Optional reference state
    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    def __init__(self, name="auxsubfieldsthermoporoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////
def auxiliary_subfields():
    """Factory associated with AuxSubfieldsThermoporoelasticity.
    """
    return AuxSubfieldsThermoporoelasticity()


# End of file
