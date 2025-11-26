# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""
Auxiliary subfields for point sources.

These subfields define the discretization for point source auxiliary fields
that may be used in visualizing or analyzing the source contribution.
"""

from pylith.topology.Subfield import Subfield


class AuxSubfieldsPointSource:
    """
    Auxiliary subfields for point moment tensor sources.
    
    This class provides subfield specifications for point source parameters
    that can be output or used in analysis.
    """

    def __init__(self, name="auxiliary_subfields"):
        """Constructor."""
        self.name = name
        self._subfields = {}

    def components(self):
        """Get list of subfield components."""
        return list(self._subfields.values())

    def addSourceLocation(self, basisOrder=0):
        """Add source location subfield specification.
        
        Args:
            basisOrder: Basis order for discretization.
        """
        subfield = Subfield("source_location")
        subfield.basisOrder = basisOrder
        self._subfields["source_location"] = subfield

    def addMomentTensor(self, basisOrder=0):
        """Add moment tensor subfield specification.
        
        Args:
            basisOrder: Basis order for discretization.
        """
        subfield = Subfield("moment_tensor")
        subfield.basisOrder = basisOrder
        self._subfields["moment_tensor"] = subfield

    def addMagnitude(self, basisOrder=0):
        """Add magnitude subfield specification.
        
        Args:
            basisOrder: Basis order for discretization.
        """
        subfield = Subfield("magnitude")
        subfield.basisOrder = basisOrder
        self._subfields["magnitude"] = subfield


# End of file
