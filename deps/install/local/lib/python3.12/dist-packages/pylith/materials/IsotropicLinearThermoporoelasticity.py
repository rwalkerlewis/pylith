# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.materials.RheologyThermoporoelasticity import RheologyThermoporoelasticity
from .materials import IsotropicLinearThermoporoelasticity as ModuleIsotropicLinearThermoporoelasticity


class IsotropicLinearThermoporoelasticity(RheologyThermoporoelasticity, ModuleIsotropicLinearThermoporoelasticity):
    """
    Isotropic linear thermoporoelastic bulk rheology.

    This rheology implements fully coupled thermoporoelasticity with:
    - Isotropic linear elasticity (drained)
    - Biot poroelasticity (pore pressure effect on stress and storage)
    - Isotropic thermal expansion (both solid and fluid)
    - Isotropic Darcy flow
    - Isotropic heat conduction

    Implements `RheologyThermoporoelasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat.bulk_rheology]
            use_reference_state = False
            
            auxiliary_subfields.biot_coefficient.basis_order = 0
            auxiliary_subfields.biot_modulus.basis_order = 0
            auxiliary_subfields.drained_bulk_modulus.basis_order = 0
            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.isotropic_permeability.basis_order = 0
            auxiliary_subfields.reference_temperature.basis_order = 0
            auxiliary_subfields.thermal_expansion_coefficient.basis_order = 0
            auxiliary_subfields.fluid_thermal_expansion.basis_order = 0
            auxiliary_subfields.thermal_conductivity.basis_order = 0
            auxiliary_subfields.specific_heat.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress and strain in computation of stress and strain."

    from .AuxSubfieldsThermoporoelasticity import AuxSubfieldsThermoporoelasticity
    auxiliarySubfields = pythia.pyre.inventory.facility("auxiliary_subfields", 
                                                         family="auxiliary_subfields", 
                                                         factory=AuxSubfieldsThermoporoelasticity)
    auxiliarySubfields.meta['tip'] = "Discretization information for physical properties and state variables."

    def __init__(self, name="isotropiclinearthermoporoelasticity"):
        """Constructor.
        """
        RheologyThermoporoelasticity.__init__(self, name)

    def preinitialize(self, problem):
        """Setup rheology.
        """
        RheologyThermoporoelasticity.preinitialize(self, problem)

        ModuleIsotropicLinearThermoporoelasticity.useReferenceState(self, self.useReferenceState)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleIsotropicLinearThermoporoelasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def thermoporoelasticity_rheology():
    """Factory associated with IsotropicLinearThermoporoelasticity.
    """
    return IsotropicLinearThermoporoelasticity()


# End of file
