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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/faults/PointSource.py
#
# @brief Python abstract base class for a point source
#
# Factory: source

from pylith.problems.Physics import Physics
from .faults import PointSource as ModulePointSource


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for fault group/nodeset/pset in mesh not specified.")
    return value


def validateLoc(value):
    """
    Validate location.
    """
    msg = "Location must be a 3 component vector (list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 3 != len(value):
        raise ValueError(msg)
    try:
        nums = list(map(float, value))
    except:
        raise ValueError(msg)
    return nums

def validateMT(value):
    """
    Validate moment tensor.
    """
    msg = " Moment tensor must be a 3D symmetric tensor (six value list)."
    if not isinstance(value, list):
        raise ValueError(msg)
    if 6 != len(value):
        raise ValueError(msg)
    try:
        nums = list(map(float, value))
    except:
        raise ValueError(msg)
    return nums

class PointSource(Physics, ModulePointSource):
    """
    Python abstract base class for a point source.

    INVENTORY

    Properties
      - *id* Source identifier
      - *label* Label identifier for source.
      - *origin_time* Origin time for point source.
      - *moment_tensor* Moment tensor of point source,
      - *dominant_frequency* Dominant frequency of Ricker function.
    Facilities
      - None

    FACTORY: fault
    """

    import pyre.inventory

    matId = pyre.inventory.int("id", default=100)
    matId.meta['tip'] = "Point source identifier (must be unique across all faults and materials)."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for point source."

    from pyre.units.time import second
    originTime = pyre.inventory.dimensional("origin_time", default=0.0 * second)
    originTime.meta['tip'] = "Origin time for earthquake rupture."

    moment_tensor = pyre.inventory.list("moment_tensor", default=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], validator=validateMT)
    moment_tensor.meta['tip'] = "Moment tensor of point source (3D symmetric tensor)."

    dominant_frequency = pyre.inventory.float("dominant_frequency", default=10.0)
    dominant_frequency.meta['tip'] = "Dominant frequency of Ricker function."


    def __init__(self, name="point_source"):
        """
        Constructor.
        """
        Physics.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Setup point source.
        """
        Physics.preinitialize(self, problem)

        ModulePointSource.setSourceId(self, self.matId)
        ModulePointSource.setSourceLabel(self, self.label)
        ModulePointSource.setMomentTensor(self, self.moment_tensor)
        ModulePointSource.setDominantFrequency(self, self.dominant_frequency)
        ModulePointSource.originTime(self, self.originTime.value)
        return

    def verifyConfiguration(self):
        """
        Verify compatibility of configuration.
        """
        return

    def _configure(self):
        """
        Setup members using inventory.
        """
        Physics._configure(self)
        return

    def _createModuleObj(self):
        """
        Create handle to corresponding C++ object.
        """
        raise NotImplementedError("Please implement _createModuleObj() in derived class.")


# End of file
