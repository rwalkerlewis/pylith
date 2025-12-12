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
# @file pylith/sources/Source.py
#
# @brief Python abstract base class for configuring sources.
#
# Factory: `source`

from pylith.problems.Physics import Physics
from .sources import Source as ModuleSource


def validateDescription(value):
    """Validate description.
    """
    if 0 == len(value):
        raise ValueError("Description for source not specified.")
    return value


class Source(Physics, ModuleSource):
    """Python abstract base class for configuring sources."""

    import pythia.pyre.inventory

    field = pythia.pyre.inventory.str("field", default="displacement")
    field.meta['tip'] = "Solution subfield associated with source."

    description = pythia.pyre.inventory.str(
        "description", default="", validator=validateDescription)
    description.meta['tip'] = "Descriptive label for source."

    labelName = pythia.pyre.inventory.str(
        "label", default="source-id",
        validator=pythia.pyre.inventory.choice(["source-id"]))
    labelName.meta['tip'] = "Name of label for source points."

    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta['tip'] = "Value of label identifying source points."

    from pylith.meshio.PointsList import PointsList
    reader = pythia.pyre.inventory.facility(
        "reader", factory=PointsList, family="points_list")
    reader.meta['tip'] = "Reader for points list."

    def __init__(self, name="source"):
        """Constructor.
        """
        Physics.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Setup source.
        """
        Physics.preinitialize(self, problem)
        ModuleSource.setSubfieldName(self, self.field)
        ModuleSource.setDescription(self, self.description)
        ModuleSource.setLabelName(self, self.labelName)
        ModuleSource.setLabelValue(self, self.labelValue)

        sourceNames, sourceCoords = self.reader.read()

        # Convert to mesh coordinate system
        from spatialdata.geocoords.Converter import convert
        convert(sourceCoords, problem.mesh().getCoordSys(), self.reader.coordsys)

        # Nondimensionalize
        if hasattr(problem.normalizer, 'lengthScale'):
            sourceCoords /= problem.normalizer.lengthScale.value
        else:
            sourceCoords /= (problem.normalizer.shearWaveSpeed.value * problem.normalizer.wavePeriod.value)

        ModuleSource.setPoints(self, sourceCoords, sourceNames)
        return


# End of file
