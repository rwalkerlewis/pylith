# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/IsotropicLinearPoroelasticity.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: poroelasticity_rheology

from .RheologyPoroelasticity import RheologyPoroelasticity
from .materials import IsotropicLinearPoroelasticity as ModuleLinearPoroelasticity


class IsotropicLinearPoroelasticity(RheologyPoroelasticity, ModuleLinearPoroelasticity):
    """
    Python material for isotropic, linearly poroelastic plane strain.

    INVENTORY

    Properties
      - *use_reference_state* Use reference stress/strain state.

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.

    FACTORY: poroelasticity_rheology
    """

    import pyre.inventory

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    from .AuxSubfieldsIsotropicLinearPoroelasticity import AuxSubfieldsIsotropicLinearPoroelasticity
    from pylith.topology.Subfield import subfieldFactory
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsIsotropicLinearPoroelasticity)
    auxiliarySubfields.meta['tip'] = "Discretization of linear poroelastic rheology physical properties."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearporoelasticity"):
        """
        Constructor.
        """
        RheologyPoroelasticity.__init__(self, name)
        return

    def preinitialize(self, mesh):
        RheologyPoroelasticity.preinitialize(self, mesh)

        ModuleLinearPoroelasticity.useReferenceState(self, self.useReferenceState)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleLinearPoroelasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def poroelasticity_rheology():
    """
    Factory associated with IsotropicLinearPoroelasticity.
    """
    return IsotropicLinearPoroelasticity()


# End of file
