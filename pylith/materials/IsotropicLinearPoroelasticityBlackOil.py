# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/IsotropicLinearPoroelasticityBlackOil.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: poroelasticityblackoil_rheology

from .RheologyMultiphasePoroelasticity import RheologyMultiphasePoroelasticity
from .materials import IsotropicLinearPoroelasticityBlackOil as ModuleLinearPoroelasticityBlackOil


class IsotropicLinearPoroelasticityBlackOil(RheologyMultiphasePoroelasticity, ModuleLinearPoroelasticityBlackOil):
    """Python material for isotropic, linearly poroelastic plane strain.

    FACTORY: multiphaseporoelasticity_rheology
    """

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useTensorPermeability = pythia.pyre.inventory.bool("use_tensor_permeability", default=False)
    useTensorPermeability.meta['tip'] = "Use tensor permeability."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearporoelasticityblackoil"):
        """Constructor.
        """
        RheologyMultiphasePoroelasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearPoroelasticityBlackOil import AuxSubfieldsIsotropicLinearPoroelasticityBlackOil
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearPoroelasticityBlackOil("auxiliary_subfields")

    def preinitialize(self, mesh):
        RheologyMultiphasePoroelasticity.preinitialize(self, mesh)

        ModuleLinearPoroelasticityBlackOil.useReferenceState(self, self.useReferenceState)
        ModuleLinearPoroelasticityBlackOil.useTensorPermeability(self, self.useTensorPermeability)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearPoroelasticityBlackOil.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def multiphaseporoelasticity_rheology():
    """Factory associated with IsotropicLinearPoroelasticityBlackOil.
    """
    return IsotropicLinearPoroelasticityBlackOil()


# End of file
