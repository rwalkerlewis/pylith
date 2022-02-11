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
# @file pylith/faults/KinSrcPoroStep.py
#
# @brief Python object for a step slip time function.
#
# Factory: eq_kinematic_src

from .KinSrcPoro import KinSrcPoro
from .faults import KinSrcPoroStep as ModuleKinSrcPoro


class KinSrcPoroStep(KinSrcPoro, ModuleKinSrcPoro):
    """Python object for a step slip time function.

    Factory: eq_kinematic_poro_src
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="kinsrcporostep"):
        """Constructor.
        """
        KinSrcPoro.__init__(self, name)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrcPoro.__init__(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_poro_src():
    """Factory associated with KinSrcPoroStep.
    """
    return KinSrcPoroStep()


# End of file
