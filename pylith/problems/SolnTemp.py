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
from .Solution import Solution as SolutionBase


class SolnTemp(PetscComponent):
    """
    Container for solution subfields with temperature subfield only.
    
    This is used for pure heat conduction problems.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnTemp
        """
    }

    import pythia.pyre.inventory

    from .SubfieldTemperature import SubfieldTemperature
    temperature = pythia.pyre.inventory.facility("temperature", family="soln_subfield", factory=SubfieldTemperature)
    temperature.meta['tip'] = "Temperature subfield."

    def __init__(self, name="solntemp"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [temperature].

        """
        return [self.temperature]


class Solution(SolutionBase):
    """Python solution field with temperature subfield.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnTemp)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
