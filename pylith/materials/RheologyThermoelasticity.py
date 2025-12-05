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


class RheologyThermoelasticity(PetscComponent):
    """
    Abstract base class for bulk rheology for thermoelasticity.
    """

    import pythia.pyre.inventory

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    auxiliaryFieldDB = pythia.pyre.inventory.facility("auxiliary_field_db", family="spatial_database", factory=SimpleDB)
    auxiliaryFieldDB.meta['tip'] = "Database for physical property parameters."

    def __init__(self, name="rheologythermoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="thermoelasticity_rheology")

    def preinitialize(self, problem):
        """Setup rheology.
        """
        self._createModuleObj()

    def addAuxiliarySubfields(self, physics, problem):
        """Add subfields for rheology to auxiliary field.
        """
        return

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        raise NotImplementedError("Implement in derived class.")


# End of file
