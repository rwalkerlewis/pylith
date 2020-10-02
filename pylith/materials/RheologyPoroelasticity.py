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
# @file pylith/materials/RheologyPoroelasticity.py
#
# @brief Python material for isotropic, linearly poroelastic, plane
# strain material.
#
# Factory: poroelasticity_rheology

from pylith.utils.PetscComponent import PetscComponent
from .materials import RheologyPoroelasticity as ModuleRheology


class RheologyPoroelasticity(PetscComponent, ModuleRheology):
    """
    Python object for bulk rheology of a poroelastic material.

    INVENTORY

    Properties
      - None

    Facilities
      - None

    FACTORY: poroelasticity_rheology
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="rheologyporoelasticity")
        return

    def preinitialize(self, mesh):
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of poroelasticity rheology '%s'." %
                           self.aliases[-1])
        self._createModuleObj()
        return

    def addAuxiliarySubfields(self, material, problem):
        for subfield in self.auxiliarySubfields.components():
            fieldName = subfield.aliases[-1]
            quadOrder = problem.defaults.quadOrder if subfield.quadOrder < 0 else subfield.quadOrder
            material.setAuxiliarySubfieldDiscretization(fieldName, subfield.basisOrder, quadOrder,
                                                        subfield.dimension, subfield.isBasisContinuous, subfield.feSpace)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Implement in derived class.")


# FACTORIES ////////////////////////////////////////////////////////////

def poroelasticity_rheology():
    """
    Factory associated with RheologyPoroelasticity.
    """
    return RheologyPoroelasticity()


# End of file
