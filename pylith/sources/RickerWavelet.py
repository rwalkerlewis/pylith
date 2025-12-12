# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/sources/RickerWavelet.py
#
# @brief Python source time functiof for a ricker wavelet.
#
# Factory: `momenttensorforce_sourcetimefunction`

from .SourceTimeFunctionMomentTensorForce import SourceTimeFunctionMomentTensorForce
from .sources import RickerWavelet as ModuleRickerWavelet


class RickerWavelet(SourceTimeFunctionMomentTensorForce, ModuleRickerWavelet):
    """Ricker wavelet source time function for moment-tensor sources."""

    import pythia.pyre.inventory

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="rickerwavelet"):
        """Constructor.
        """
        SourceTimeFunctionMomentTensorForce.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsSourceTime import AuxSubfieldsSourceTime
        self.auxiliarySubfields = AuxSubfieldsSourceTime("auxiliary_subfields")

    def preinitialize(self, problem):
        SourceTimeFunctionMomentTensorForce.preinitialize(self, problem)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleRickerWavelet.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def momenttensorforce_sourcetimefunction():
    """Factory associated with RickerWavelet.
    """
    return RickerWavelet()


# End of file
