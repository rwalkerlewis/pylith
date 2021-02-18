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

#include "pylith/faults/PointSource.hh" // implementation of object methods

#include "pylith/faults/TopologyOps.hh" // USES TopologyOps
#include "pylith/materials/DerivedFactoryElasticity.hh" // USES DerivedFactoryElasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::checkTopology()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/PointSource.hh" // USES PointSource

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
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;

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
              * @param[in] formulation Formulation for equations.              
              */
              static
              void setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                         const pylith::faults::PointSource& pointSource,
                                         const pylith::topology::Field& solution,
                                         const pylith::problems::Physics::FormulationEnum formulation);

              static const char* pyreComponent;

           }; // _PointSource
           const char* _PointSource::pyreComponent = "pointsource";

       } // faults
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default Constructor
pylith::faults::PointSource::PointSource(void) :
    _sourceId(100),
    _sourceLabel(""),
    _momentTensorConvention("harvard"),
    _originTime(0.0) { // t
    _momentTensor[0] = 0.0; // Mrr / Mxx
    _momentTensor[1] = 0.0; // Mtt / Myy
    _momentTensor[2] = 0.0; // Mpp / Mzz
    _momentTensor[3] = 0.0; // Mrt / Mxy
    _momentTensor[4] = 0.0; // Mrp / Mxz
    _momentTensor[5] = 0.0; // Mtp / Myz
    _location[0] = 0.0; // x
    _location[1] = 0.0; // y
    _location[2] = 0.0; // z
    _dominantFrequency = 10; // f_0, Hz
    pylith::utils::PyreComponent::setName(_PointSource::pyreComponent);
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
// Set identifier for point source cells.
void
pylith::faults::PointSource::setSourceId(const int value) {
    PYLITH_COMPONENT_DEBUG("setSourceId(value="<<value<<")");

    _sourceId = value;
} // setPointSourceId


// ---------------------------------------------------------------------------------------------------------------------
// Get identifier for source cells.
int
pylith::faults::PointSource::getSourceId(void) const {
    return _sourceId;
} // getPointSourceId

// ---------------------------------------------------------------------------------------------------------------------
// Set label marking point source
void
pylith::faults::PointSource::setSourceLabel(const char* value) {
    PYLITH_COMPONENT_DEBUG("setSourceLabel(value="<<value<<")");

    _sourceLabel = value;
} // setSurfaceMarkerLabel

// ---------------------------------------------------------------------------------------------------------------------
// Get label marking point source.
const char*
pylith::faults::PointSource::getSourceLabel(void) const {
    return _sourceLabel.c_str();
} // getPointSourceLabel

// ---------------------------------------------------------------------------------------------------------------------
// Set origin time for point source
void
pylith::faults::PointSource::setOriginTime(const PylithReal value) {
    PYLITH_COMPONENT_DEBUG("setPointShift(t="<<value<<")");

    _originTime = value;

} // setOriginTime

// ----------------------------------------------------------------------
// Get origin time for point source.
PylithReal
pylith::faults::PointSource::getOriginTime(void) const {
    return _originTime;
} // getOriginTime

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
// Set dominant frequency of Ricker function
void
pylith::faults::PointSource::setDominantFrequency(const PylithReal value) {
    PYLITH_COMPONENT_DEBUG("setDominantFrequency(f="<<value<<")");

    _dominantFrequency = value;

} // setDominantFrequency

// ---------------------------------------------------------------------------------------------------------------------
// Set dominant frequency of Ricker function
PylithReal
pylith::faults::PointSource::getDominantFrequency(const PylithReal value) {
    return _dominantFrequency;
} // getDominantFrequency

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
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::faults::PointSource::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setMaterialId(getSourceId());
    //integrator->setSourceLabel(getSourceLabel());

    _PointSource::setKernelsRHSResidual(integrator, *this, solution, _formulation);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator

// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::faults::PointSource::createDerivedField(const pylith::topology::Field& solution,
                                                     const pylith::topology::Mesh& domainMesh) {
    return NULL;
} // createDerivedField

// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::faults::PointSource::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelConstants(dt="<<dt<<")");

    if (11 != _kernelConstants.size()) { _kernelConstants.resize(11);}
    _kernelConstants[0] = _momentTensor[0]; // Mrr / Mxx
    _kernelConstants[1] = _momentTensor[1]; // Mtt / Myy
    _kernelConstants[2] = _momentTensor[2]; // Mpp / Mzz
    _kernelConstants[3] = _momentTensor[3]; // Mrt / Mxy
    _kernelConstants[4] = _momentTensor[4]; // Mrp / Mxz
    _kernelConstants[5] = _momentTensor[5]; // Mtp / Myz
    _kernelConstants[6] = _location[0];
    _kernelConstants[7] = _location[1];
    _kernelConstants[8] = _location[2];
    _kernelConstants[9] = _originTime;
    _kernelConstants[10] = _dominantFrequency;

    PYLITH_METHOD_END;
} // _updateKernelConstants

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for RHS residual.
void
pylith::faults::_PointSource::setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                         const pylith::faults::PointSource& pointSource,
                                                         const pylith::topology::Field& solution,
                                                         const pylith::problems::Physics::FormulationEnum formulation) {

    PYLITH_METHOD_BEGIN;
    journal::debug_t debug(_PointSource::pyreComponent);
    debug << journal::at(__HERE__)
          << "setKernelsRHSResidual(integrator="<<integrator<<", pointsource="<<typeid(pointSource).name()<<", solution="
          << solution.getLabel()<<")"
          << journal::endl;

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();

    std::vector<ResidualKernels> kernels;


    PetscPointFunc g0u = NULL;
    PetscPointFunc g1u = NULL;
    PetscPointFunc g0v = NULL;
    PetscPointFunc g1v = NULL;


    // Assign specific kernel based on source / time function selected

    g0u = pylith::fekernels::PointSource::g0u_ricker;
    kernels.resize(2);
    kernels[0] = ResidualKernels("displacement", g0u, g1u);
    kernels[1] = ResidualKernels("velocity", g0v, g1v);


    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // setKernelsRHSResidual


// End of file
