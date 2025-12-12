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


class SolnPresTemp(PetscComponent):
    """
    Container for solution subfields with pressure and temperature subfields.
    
    This is used for coupled fluid flow and heat transfer problems.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnPresTemp
        """
    }

    import pythia.pyre.inventory

    from .SubfieldPressure import SubfieldPressure
    pressure = pythia.pyre.inventory.facility("pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."

    from .SubfieldTemperature import SubfieldTemperature
    temperature = pythia.pyre.inventory.facility("temperature", family="soln_subfield", factory=SubfieldTemperature)
    temperature.meta['tip'] = "Temperature subfield."

    def __init__(self, name="solnprestemp"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [pressure, temperature].

        """
        return [self.pressure, self.temperature]


class Solution(SolutionBase):
    """Python solution field with pressure and temperature subfields.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnPresTemp)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
