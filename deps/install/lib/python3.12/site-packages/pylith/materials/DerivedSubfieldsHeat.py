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


class DerivedSubfieldsHeat(PetscComponent):
    """
    Derived subfields associated with the heat equation.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_heat.derived_fields]
            heat_flux.basis_order = 0
        """
    }
    """Python container for heat equation derived subfields.

    FACTORY: derived_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    heatFlux = pythia.pyre.inventory.facility("heat_flux", family="derived_subfield", factory=Subfield)
    heatFlux.meta['tip'] = "Heat flux subfield."

    def __init__(self, name="derivedsubfieldsheat"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="derived_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def derived_subfields():
    """Factory associated with DerivedSubfieldsHeat.
    """
    return DerivedSubfieldsHeat()


# End of file
