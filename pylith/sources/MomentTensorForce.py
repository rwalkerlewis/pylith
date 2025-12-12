# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
#
# @file pylith/sources/MomentTensorForce.py
#
# @brief Python component for a moment-tensor point source.
#
# Factory: `source`

from pylith.sources.TimeHistoryWavelet import TimeHistoryWavelet
from .Source import Source
from .sources import MomentTensorForce as ModuleMomentTensorForce


class MomentTensorForce(Source, ModuleMomentTensorForce):
    """Moment-tensor point source."""

    import pythia.pyre.inventory

    source_time_function = pythia.pyre.inventory.facility(
        "source_time_function", family="momenttensorforce_sourcetimefunction", factory=TimeHistoryWavelet)
    source_time_function.meta['tip'] = "Source time function for moment tensor force."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="momenttensorforce"):
        """Constructor.
        """
        Source.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsMomentTensorForce import AuxSubfieldsMomentTensorForce
        self.auxiliarySubfields = AuxSubfieldsMomentTensorForce(
            "auxiliary_subfields")

    def preinitialize(self, problem):
        """Setup source.
        """
        self.source_time_function.preinitialize(problem)
        Source.preinitialize(self, problem)

        self.source_time_function.addAuxiliarySubfields(self, problem)

        return

    def _createModuleObj(self):
        """Create handle to C++ MomentTensorForce.
        """
        ModuleMomentTensorForce.__init__(self)
        # Material sets auxiliary db in source_time_function.
        ModuleMomentTensorForce.setSourceTimeFunction(
            self, self.source_time_function)
        return


# Factories

def source():
    """Factory associated with MomentTensorForce.
    """
    return MomentTensorForce()


# End of file
