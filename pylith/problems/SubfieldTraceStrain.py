#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/problems/SubfieldTraceStrain.py
##
# @brief Python object for trace_strain subfield.
##
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield


class SubfieldTraceStrain(SolutionSubfield):
    """
    Python object for trace_strain subfield.

    INVENTORY

    Properties
      - *alias* User-specified name for subfield.

    Facilities
      - None

    FACTORY: subfield
    """

    import pyre.inventory

    from .SolutionSubfield import validateAlias
    userAlias = pyre.inventory.str("alias", default="trace_strain", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    fieldName = "trace_strain"

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfieldtracestrain"):
        """
        Constructor.
        """
        SolutionSubfield.__init__(self, name)
        return

    def initialize(self, normalizer, spaceDim):
        """
        Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        from pyre.units.unit import one
        self.vectorFieldType = Field.SCALAR
        self.scale = one
        self._setComponents(spaceDim)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        SolutionSubfield._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """
    Factory associated with SubfieldTraceStrain.
    """
    return SubfieldTraceStrain()


# End of file
