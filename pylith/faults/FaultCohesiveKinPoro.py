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

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveKinPoro as ModuleFaultCohesiveKinPoro


def eqsrcFactory(name):
    """
    Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .KinSrcPoroStep import KinSrcPoroStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcPoroStep)


class FaultCohesiveKinPoro(FaultCohesive, ModuleFaultCohesiveKinPoro):
    """Python object for a fault surface with kinematic (prescribed) slip and fluid diffusion
    implemented with cohesive elements.

    FACTORY: fault
    """
    DOC_CONFIG = {
        "cfg": """
            # Specify prescribed slip on a fault via two earthquakes in a 2D domain.
            [pylithapp.problem.interfaces.fault]
            label = fault
            edge = fault_edge

            observers.observer.data_fields = [slip]

            # Two earthquakes with different slip time functions.
            eq_ruptures = [quake10, quake50]
            quake10 = pylith.faults.KinSrcPoroBrune
            quake50 = pylith.faults.KinSrcPoroLiuCosine

            # Rupture parameters for the first earthquake.
            [pylithapp.problem.interfaces.fault.eq_ruptures.quake10]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Fault rupture auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

            # Rupture parameters for the second earthquake.
            [pylithapp.problem.interfaces.fault.eq_ruptures.quake50]
            origin_time = 50*year
            
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Fault rupture auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -1.0*m, 0.0*m]
            """
    }

    import pythia.pyre.inventory

    from .SingleRupturePoro import SingleRupturePoro
    eqRuptures = pythia.pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupturePoro)
    eqRuptures.meta['tip'] = "Kinematic poroelastic earthquake sources information."

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)
    auxiliaryFieldDB.meta['tip'] = "Database for fault thickness, porosity, beta_p, beta_sigma, permeability, fluid viscosity, and slip."

    def __init__(self, name="faultcohesivekinporo"):
        """Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsFaultPoro import AuxSubfieldsFaultPoro
        self.auxiliarySubfields = AuxSubfieldsFaultPoro("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_is_root
        if mpi_is_root():
            self._info.log("Pre-initializing fault '%s'." % self.labelName)

        FaultCohesive.preinitialize(self, problem)

        for eqsrc in self.eqRuptures.components():
            eqsrc.preinitialize()
        ModuleFaultCohesiveKinPoro.setEqRuptures(
            self, self.eqRuptures.inventory.facilityNames(), self.eqRuptures.components())

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