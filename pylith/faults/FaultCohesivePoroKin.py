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

from .FaultCohesivePoro import FaultCohesivePoro
from .faults import FaultCohesivePoroKin as ModuleFaultCohesivePoroKin


def eqsrcFactory(name):
    """
    Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .KinSrcPoroStep import KinSrcPoroStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcPoroStep)


class FaultCohesivePoroKin(FaultCohesivePoro, ModuleFaultCohesivePoroKin):
    """Python object for a fault surface with kinematic (prescribed) slip and fluid diffusion
    implemented with cohesive elements.

    FACTORY: fault
    """
    DOC_CONFIG = {
        """
        @file pylith/faults/FaultCohesivePoroKin.py
        
        @brief Python object for a poroelastic fault surface with kinematic
        (prescribed) slip, fluid diffusion implemented with cohesive elements.
        
        Factory: fault
        """
    }  

    import pythia.pyre.inventory

    from .SingleRupturePoro import SingleRupturePoro
    eqRuptures = pythia.pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupturePoro)
    eqRuptures.meta['tip'] = "Kinematic poroelastic earthquake sources information."

    from pylith.utils.NullComponent import NullComponent
    auxFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)
    auxFieldDB.meta['tip'] = "Database for fault thickness, porosity, beta_p, beta_sigma, permeability, fluid viscosity, and slip."

    def __init__(self, name="faultcohesiveporokin"):
        """Initialize configuration.
        """
        FaultCohesivePoro.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Pre-initializing fault '%s'." % self.labelName)

        FaultCohesivePoro.preinitialize(self, problem)

        for eqsrc in self.eqRuptures.components():
            eqsrc.preinitialize()
        ModuleFaultCohesivePoroKin.setEqRuptures(
            self, self.eqRuptures.inventory.facilityNames(), self.eqRuptures.components())

        ModuleFaultCohesivePoroKin.auxFieldDB(self, self.auxFieldDB)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesivePoro.verifyConfiguration(self)
        ModuleFaultCohesivePoroKin.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqRuptures.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """Cleanup.
        """
        for eqsrc in self.eqRuptures.components():
            eqsrc.finalize()
        FaultCohesivePoro.finalize(self)
        return

    def _configure(self):
        """Setup members using inventory.
        """
        FaultCohesivePoro._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesivePoroKin.
        """
        ModuleFaultCohesivePoroKin.__init__(self)
        return


# Factories

def fault():
    """Factory associated with FaultCohesivePoroKin.
    """
    return FaultCohesivePoroKin()


# End of file