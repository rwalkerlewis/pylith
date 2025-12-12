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


class RheologyHeat(PetscComponent):
    """
    Abstract base class for bulk rheology for heat conduction.
    """

    import pythia.pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxiliaryFieldDB = pythia.pyre.inventory.facility(
        "db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    auxiliaryFieldDB.meta['tip'] = "Database for physical property parameters."

    def __init__(self, name):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="heat_rheology")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        self._createModuleObj()

    def addAuxiliarySubfields(self, material, problem):
        """Add subfields to auxiliary field.
        """
        pass

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        raise NotImplementedError("Implement in derived class.")


# End of file
