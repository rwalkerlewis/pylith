# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .RheologyThermoelasticity import RheologyThermoelasticity
from .materials import IsotropicLinearThermoelasticity as ModuleIsotropicLinearThermoelasticity


class IsotropicLinearThermoelasticity(RheologyThermoelasticity, ModuleIsotropicLinearThermoelasticity):
    """
    Isotropic, linear thermoelastic bulk rheology.

    This rheology provides:
    - Linear elastic stress-strain relationship with thermal strain
    - Isotropic thermal conductivity (Fourier's law)
    - Thermal coupling: stress = C:(strain - alpha*(T-Tref)*I)

    Implements `RheologyThermoelasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_thermoelastic.bulk_rheology]
            use_reference_state = False

            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.bulk_modulus.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    def __init__(self, name="isotropiclinearthermoelasticity"):
        """Constructor.
        """
        RheologyThermoelasticity.__init__(self, name)

    def preinitialize(self, problem):
        """Setup rheology.
        """
        RheologyThermoelasticity.preinitialize(self, problem)
        ModuleIsotropicLinearThermoelasticity.useReferenceState(self, self.useReferenceState)

    def addAuxiliarySubfields(self, physics, problem):
        """Add subfields for rheology to auxiliary field.
        """
        from pylith.topology.Field import Field
        from .AuxSubfieldsThermoelasticity import AuxSubfieldsThermoelasticity

        auxiliarySubfields = AuxSubfieldsThermoelasticity("auxiliary_subfields")
        self._setupAuxiliarySubfields(physics, auxiliarySubfields, problem)

    def _setupAuxiliarySubfields(self, physics, auxiliarySubfields, problem):
        """Configure auxiliary subfields.
        """
        # Add shear modulus and bulk modulus for elasticity
        # These are typically queried via vs/vp from the database
        for subfieldName in ["shear_modulus", "bulk_modulus"]:
            if hasattr(auxiliarySubfields, subfieldName):
                subfield = getattr(auxiliarySubfields, subfieldName)
                physics.setAuxiliarySubfieldDiscretization(
                    subfieldName,
                    subfield.basisOrder,
                    subfield.quadOrder,
                    problem.defaults.quadrature_order,
                    subfield.isBasisContinuous,
                    subfield.feSpace
                )

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleIsotropicLinearThermoelasticity.__init__(self)


# Factories

def thermoelasticity_rheology():
    """Factory associated with IsotropicLinearThermoelasticity.
    """
    return IsotropicLinearThermoelasticity()


# End of file
