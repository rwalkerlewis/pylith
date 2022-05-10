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
# @file pylith/faults/FaultCohesiveKinPoro.py
#
# @brief Python object for a poroelastic fault surface with kinematic
# (prescribed) slip, fluid diffusion implemented with cohesive elements.
#
# Factory: fault

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveKinPoro as ModuleFaultCohesiveKinPoro

# ITEM FACTORIES ///////////////////////////////////////////////////////


def eqsrcFactory(name):
    """Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .KinSrcPoroStep import KinSrcPoroStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcPoroStep)


class FaultCohesiveKinPoro(FaultCohesive, ModuleFaultCohesiveKinPoro):
    """Python object for a fault surface with kinematic (prescribed) slip and fluid diffusion
    implemented with cohesive elements.

    FACTORY: fault
    """

    import pythia.pyre.inventory

    from .SingleRupturePoro import SingleRupturePoro
    eqRuptures = pythia.pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupturePoro)
    eqRuptures.meta['tip'] = "Kinematic earthquake sources information."

    from pylith.utils.NullComponent import NullComponent
    auxFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)
    auxFieldDB.meta['tip'] = "Database for fault thickness, porosity, beta_p, beta_sigma, permeability, fluid viscosity, and slip."

    def __init__(self, name="faultcohesivekinporo"):
        """Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Pre-initializing fault '%s'." % self.labelName)

        FaultCohesive.preinitialize(self, problem)

        for eqsrc in self.eqRuptures.components():
            eqsrc.preinitialize()
        print(self.eqRuptures.components())
        ModuleFaultCohesiveKinPoro.setEqRuptures(
            self, self.eqRuptures.inventory.facilityNames(), self.eqRuptures.components())

        ModuleFaultCohesiveKinPoro.auxFieldDB(self, self.auxFieldDB)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveKinPoro.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqRuptures.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """Cleanup.
        """
        for eqsrc in self.eqRuptures.components():
            eqsrc.finalize()
        FaultCohesive.finalize(self)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        FaultCohesive._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesiveKinPoro.
        """
        ModuleFaultCohesiveKinPoro.__init__(self)
        return


# Factories

def fault():
    """Factory associated with FaultCohesiveKinPoro.
    """
    return FaultCohesiveKinPoro()


# End of file