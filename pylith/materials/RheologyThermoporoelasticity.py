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


class RheologyThermoporoelasticity(PetscComponent):
    """
    Abstract base class for bulk rheology for thermoporoelasticity.

    This defines the interface for thermoporoelastic rheologies that combine:
    - Elastic behavior (stress-strain relationship)
    - Pore fluid effects (Biot poroelasticity)
    - Thermal effects (thermal expansion, heat conduction)
    """

    import pythia.pyre.inventory

    def __init__(self, name="rheologythermoporoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="thermoporoelasticity_rheology")

    def preinitialize(self, problem):
        """Setup rheology.
        """
        self._createModuleObj()

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        raise NotImplementedError("Implement in child class.")


# End of file
