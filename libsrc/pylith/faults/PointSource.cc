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
              void setKernelsLHSResidual(pylith::feassemble::IntegratorBoundary* integrator,
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
    _momentTensor[0] = 0.0; // Mrr
    _momentTensor[1] = 0.0; // Mtt
    _momentTensor[2] = 0.0; // Mpp
    _momentTensor[3] = 0.0; // Mrt
    _momentTensor[4] = 0.0; // Mrp
    _momentTensor[5] = 0.0; // Mtp

    _pointLocation[0] = 0.0; // x
    _pointLocation[1] = 0.0; // y
    _pointLocation[2] = 0.0; // z

    _pointElapsed = 0.0; // t
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


// End of file
