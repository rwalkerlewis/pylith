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


class SolnDispPresTracStrainTempLagrange(PetscComponent):
    """Container for solution subfields for thermoporoelasticity with faults.

    Includes displacement, pressure, trace_strain, temperature, and fault Lagrange multiplier.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnDispPresTracStrainTempLagrange
        """
    }

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement

    displacement = pythia.pyre.inventory.facility(
        "displacement", family="soln_subfield", factory=SubfieldDisplacement
    )
    displacement.meta["tip"] = "Displacement subfield."

    from .SubfieldPressure import SubfieldPressure

    pressure = pythia.pyre.inventory.facility(
        "pressure", family="soln_subfield", factory=SubfieldPressure
    )
    pressure.meta["tip"] = "Pressure subfield."

    from .SubfieldTraceStrain import SubfieldTraceStrain

    traceStrain = pythia.pyre.inventory.facility(
        "trace_strain", family="soln_subfield", factory=SubfieldTraceStrain
    )
    traceStrain.meta["tip"] = "Trace strain subfield."

    from .SubfieldTemperature import SubfieldTemperature

    temperature = pythia.pyre.inventory.facility(
        "temperature", family="soln_subfield", factory=SubfieldTemperature
    )
    temperature.meta["tip"] = "Temperature subfield."

    from .SubfieldLagrangeFault import SubfieldLagrangeFault

    lagrangeFault = pythia.pyre.inventory.facility(
        "lagrange_multiplier_fault",
        family="soln_subfield",
        factory=SubfieldLagrangeFault,
    )
    lagrangeFault.meta["tip"] = "Fault Lagrange multiplier subfield."

    def __init__(self, name="SolnDispPresTracStrainTempLagrange"):
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Return components in a deterministic order."""
        return [
            self.displacement,
            self.pressure,
            self.traceStrain,
            self.temperature,
            self.lagrangeFault,
        ]


class Solution(SolutionBase):
    """Python solution field for thermoporoelasticity with faults."""

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory

    subfields = pythia.pyre.inventory.facilityArray(
        "subfields",
        family="soln_subfields",
        itemFactory=subfieldFactory,
        factory=SolnDispPresTracStrainTempLagrange,
    )
    subfields.meta["tip"] = "Subfields in solution."


def solution():
    """Factory associated with Solution."""

    return Solution()


# End of file
