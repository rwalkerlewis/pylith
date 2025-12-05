# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .RheologyHeat import RheologyHeat
from .materials import IsotropicHeat as ModuleIsotropicHeat


class IsotropicHeat(RheologyHeat, ModuleIsotropicHeat):
    """
    Isotropic heat conduction rheology.

    This rheology implements Fourier's law for isotropic heat conduction:
    q = -k ∇T

    where:
    - q is the heat flux vector
    - k is the thermal conductivity (scalar for isotropic)
    - T is the temperature

    Implements `RheologyHeat`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_heat.bulk_rheology]
            db_auxiliary_field = spatialdata.spatialdb.SimpleDB
            db_auxiliary_field.description = Heat material properties
            db_auxiliary_field.iohandler.filename = mat_heat.spatialdb
        """
    }

    def __init__(self, name="isotropicheat"):
        """Constructor.
        """
        RheologyHeat.__init__(self, name)

    def _createModuleObj(self):
        """Create handle to C++ IsotropicHeat.
        """
        ModuleIsotropicHeat.__init__(self)


# Factories

def heat_rheology():
    """Factory associated with IsotropicHeat.
    """
    return IsotropicHeat()


# End of file
