// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "PointSource.hh" // implementation of object methods

#include "pylith/faults/TopologyOps.hh" // USES TopologyOps

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/MeshOps.hh" // USES MeshOps::checkTopology()

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensionalizer
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _PointSource {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Set kernels for LHS residual.
              *
              * @param[out] integrator Integrator for point
              * @param[in] pointSource PointSource for point source condition.
              * @param[in] solution Solution field.
              */
              static
              void setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                         const pylith:faults:PointSource& pointSource,
                                         const pylith::topology::Field& solution);
              static const char* pyreComponent;

           }; // _PointSource
           const char* _PointSource::pyreComponent = "pointsource";

       } // faults
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default Constructor
pylith::faults::PointSource::PointSource(void) :
    _momentTensorConvention(""),
    _sourceLabel("") {
    _momentTensor[0] = 0.0; // Mrr / Mxx
    _momentTensor[1] = 0.0; // Mtt / Myy
    _momentTensor[2] = 0.0; // Mpp / Mzz
    _momentTensor[3] = 0.0; // Mrt / Mxy
    _momentTensor[4] = 0.0; // Mrp / Mxz
    _momentTensor[5] = 0.0; // Mtp / Myz

    _pointLocation[0] = 0.0; // x
    _pointLocation[1] = 0.0; // y
    _pointLocation[2] = 0.0; // z

    _pointShift = 0.0; // t

    _dominantFrequency = 10; // f_0, Hz
} // constructor

// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::faults::PointSource::~PointSource(void) {
  deallocate();
} // destructor

// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures
void
pylith::faults::PointSource::deallocate(void) {
  PYLITH_METHOD_BEGIN;

  pylith::problems::Physics::deallocate();

  PYLITH_METHOD_END;
} // deallocate

// ---------------------------------------------------------------------------------------------------------------------
// Set moment tensor parameters
void
pylith::faults::PointSource::setMomentTensor(const PylithReal vec[6]) {
    // Given Harvard CMT convention
    PYLITH_COMPONENT_DEBUG("setMomentTensor(Mrr="<<vec[0]<<", Mtt="<<vec[1]<<", Mpp="<<vec[2]
                                       <<", Mrt="<<vec[3]<<", Mrp="<<vec[4]<<", Mtp="<<vec[5]<<")");

     for (int i = 0; i < 6; ++i) {
        _momentTensor[i] = vec[i];
     } // for
} // setMomentTensor

// ---------------------------------------------------------------------------------------------------------------------
// Set location of point source
void
pylith::faults::PointSource::setPointLocation(const PylithReal vec[3]) {
    PYLITH_COMPONENT_DEBUG("setPointLocation(x="<<vec[0]<<", y="<<vec[1]<<", z="<<vec[2]<<")");

     for (int i = 0; i < 3; ++i) {
        _pointLocation[i] = vec[i];
     } // for
} // setPointLocation

// ---------------------------------------------------------------------------------------------------------------------
// Set time delay of point source
void
pylith::faults::PointSource::setPointShift(const PylithReal vec[1]]) {
    PYLITH_COMPONENT_DEBUG("setPointShift(t="<<vec[0]<<")");

    _pointShift = vec;

} // setPointShift

// ---------------------------------------------------------------------------------------------------------------------
// Set dominant frequency of Ricker function
void
pylith::faults::PointSource::setDominantFrequency(const PylithReal vec[1]]) {
    PYLITH_COMPONENT_DEBUG("setDominantFrequency(f="<<vec[0]<<")");

    _dominantFrequency = vec;

} // setDominantFrequency

// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::PointSource::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'Elasticity'.");
    } // if
    switch (_formulation) {
    case QUASISTATIC:
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity")) {
            throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'Elasticity' using "
                                     "'dynamic' or 'dynamic_imex' formulation.");
        } // if
        break;
    default:
        PYLITH_COMPONENT_FIREWALL("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration

// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::PointSource::_getAuxiliaryFactory(void) {
    return NULL;
} // _getAuxiliaryFactory

// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::PointSource::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt="<<dt<<")");

    if (10 != _kernelConstants.size()) { _kernelConstants.resize(11);}
    _kernelConstants[0] = _momentTensor[0]; // Mrr / Mxx
    _kernelConstants[1] = _momentTensor[1]; // Mtt / Myy
    _kernelConstants[2] = _momentTensor[2]; // Mpp / Mzz
    _kernelConstants[3] = _momentTensor[3]; // Mrt / Mxy
    _kernelConstants[4] = _momentTensor[4]; // Mrp / Mxz
    _kernelConstants[5] = _momentTensor[5]; // Mtp / Myz
    _kernelConstants[6] = _pointLocation[0];
    _kernelConstants[7] = _pointLocation[1];
    _kernelConstants[8] = _pointLocation[2];
    _kernelConstants[9] = _pointShift;
    _kernelConstants[10] = _dominantFrequency;

    PYLITH_METHOD_END;
} // _updateKernelConstants

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual.
void
pylith::faults::_FaultCohesiveKin::setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                         const pylith::faults::PointSource& pointSource,
                                                         const pylith::topology::Field& solution) {
                                                             
    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_PointSource::pyreComponent);
    debug << journal::at(__HERE__)
          << "setKernelsRHSResidual(integrator="<<integrator<<", pointsource="<<typeid(pointSource).name()<<", solution="
          << solution.getLabel()<<")"
          << journal::endl;

    PetscPointFunc g0 = NULL;
    PetscPointFunc g1 = NULL;

    // Assign specific kernel based on source / time function selected

    g0 = pylith::fekernels::PointSource::g0u_ricker;

    std::vector<ResidualKernels> kernels(1);
    kernels[0] = ResidualKernels(pointSource.getSubfieldName(), g0, g1);

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END
} // setKernelsRHSResidual


// End of file
