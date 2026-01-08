# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
# @file pylith/materials/IsotropicLinearBlackOilPoroelasticity.py
#
# @brief Python material for isotropic, linear poroelasticity with
# black oil formulation for pressure-dependent fluid properties.
#
# Factory: poroelasticity_rheology

from .RheologyPoroelasticity import RheologyPoroelasticity
from .materials import IsotropicLinearBlackOilPoroelasticity as ModuleBlackOilPoroelasticity


class IsotropicLinearBlackOilPoroelasticity(RheologyPoroelasticity, ModuleBlackOilPoroelasticity):
    """
    Isotropic linear poroelasticity with black oil formulation.

    The black oil model extends standard poroelasticity with pressure-dependent
    fluid properties commonly used in petroleum reservoir simulation:
    - Fluid compressibility varies with pressure
    - Fluid viscosity varies with pressure

    Implements `RheologyPoroelasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_blackoil.rheology]
            use_reference_state = False
            use_tensor_permeability = False

            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.drained_bulk_modulus.basis_order = 0
            auxiliary_subfields.reference_pressure.basis_order = 0
            auxiliary_subfields.fluid_compressibility.basis_order = 0
            auxiliary_subfields.fluid_compressibility_coefficient.basis_order = 0
            auxiliary_subfields.viscosity_coefficient.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useTensorPermeability = pythia.pyre.inventory.bool("use_tensor_permeability", default=False)
    useTensorPermeability.meta['tip'] = "Use tensor permeability."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearblackoilporoelasticity"):
        """Constructor.
        """
        RheologyPoroelasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearBlackOilPoroelasticity import AuxSubfieldsIsotropicLinearBlackOilPoroelasticity
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearBlackOilPoroelasticity("auxiliary_subfields")

        from .DerivedSubfieldsPoroelasticity import DerivedSubfieldsPoroelasticity
        self.derivedSubfields = DerivedSubfieldsPoroelasticity("derived_subfields")

    def preinitialize(self, mesh):
        RheologyPoroelasticity.preinitialize(self, mesh)

        ModuleBlackOilPoroelasticity.useReferenceState(self, self.useReferenceState)
        ModuleBlackOilPoroelasticity.useTensorPermeability(self, self.useTensorPermeability)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleBlackOilPoroelasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def poroelasticity_rheology():
    """Factory associated with IsotropicLinearBlackOilPoroelasticity.
    """
    return IsotropicLinearBlackOilPoroelasticity()


# End of file
