// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/sources/PointForce.hh" // implementation of object methods

#include "pylith/sources/PointForceAuxiliaryFactory.hh" // USES PointForceAuxiliaryFactory
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/PointForce.hh" // USES PointForce kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::PointForce::PointForce(void) :
    _magnitude(1.0e+15),
    _originTime(0.0),
    _dominantFrequency(1.0),
    _timeDelay(1.0),
    _auxiliaryFactory(new pylith::sources::PointForceAuxiliaryFactory) {
    pylith::utils::PyreComponent::setName("pointforce");

    // Default isotropic moment tensor
    _momentTensor.resize(6);
    _momentTensor[0] = 1.0; // Mxx
    _momentTensor[1] = 1.0; // Myy
    _momentTensor[2] = 1.0; // Mzz
    _momentTensor[3] = 0.0; // Mxy
    _momentTensor[4] = 0.0; // Mxz
    _momentTensor[5] = 0.0; // Myz

    // Default location at origin
    _location.resize(3);
    _location[0] = 0.0;
    _location[1] = 0.0;
    _location[2] = 0.0;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::PointForce::~PointForce(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::sources::PointForce::deallocate(void) {
    Physics::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set source location.
void
pylith::sources::PointForce::setLocation(const PylithReal* location,
                                         const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setLocation(location="<<location<<", size="<<size<<")");

    assert(location);
    assert(size > 0 && size <= 3);

    _location.resize(size);
    for (int i = 0; i < size; ++i) {
        _location[i] = location[i];
    }

    PYLITH_METHOD_END;
} // setLocation


// ------------------------------------------------------------------------------------------------
// Get source location.
const PylithReal*
pylith::sources::PointForce::getLocation(void) const {
    return _location.data();
} // getLocation


// ------------------------------------------------------------------------------------------------
// Get spatial dimension of source location.
int
pylith::sources::PointForce::getLocationSize(void) const {
    return _location.size();
} // getLocationSize


// ------------------------------------------------------------------------------------------------
// Set moment tensor components.
void
pylith::sources::PointForce::setMomentTensor(const PylithReal* momentTensor,
                                             const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setMomentTensor(momentTensor="<<momentTensor<<", size="<<size<<")");

    assert(momentTensor);
    assert(size == 3 || size == 4 || size == 6);

    _momentTensor.resize(size);
    for (int i = 0; i < size; ++i) {
        _momentTensor[i] = momentTensor[i];
    }

    PYLITH_METHOD_END;
} // setMomentTensor


// ------------------------------------------------------------------------------------------------
// Get moment tensor components.
const PylithReal*
pylith::sources::PointForce::getMomentTensor(void) const {
    return _momentTensor.data();
} // getMomentTensor


// ------------------------------------------------------------------------------------------------
// Get number of moment tensor components.
int
pylith::sources::PointForce::getMomentTensorSize(void) const {
    return _momentTensor.size();
} // getMomentTensorSize


// ------------------------------------------------------------------------------------------------
// Set source magnitude.
void
pylith::sources::PointForce::setMagnitude(const PylithReal value) {
    _magnitude = value;
} // setMagnitude


// ------------------------------------------------------------------------------------------------
// Get source magnitude.
PylithReal
pylith::sources::PointForce::getMagnitude(void) const {
    return _magnitude;
} // getMagnitude


// ------------------------------------------------------------------------------------------------
// Set origin time.
void
pylith::sources::PointForce::setOriginTime(const PylithReal value) {
    _originTime = value;
} // setOriginTime


// ------------------------------------------------------------------------------------------------
// Get origin time.
PylithReal
pylith::sources::PointForce::getOriginTime(void) const {
    return _originTime;
} // getOriginTime


// ------------------------------------------------------------------------------------------------
// Set dominant frequency.
void
pylith::sources::PointForce::setDominantFrequency(const PylithReal value) {
    _dominantFrequency = value;
} // setDominantFrequency


// ------------------------------------------------------------------------------------------------
// Get dominant frequency.
PylithReal
pylith::sources::PointForce::getDominantFrequency(void) const {
    return _dominantFrequency;
} // getDominantFrequency


// ------------------------------------------------------------------------------------------------
// Set time delay.
void
pylith::sources::PointForce::setTimeDelay(const PylithReal value) {
    _timeDelay = value;
} // setTimeDelay


// ------------------------------------------------------------------------------------------------
// Get time delay.
PylithReal
pylith::sources::PointForce::getTimeDelay(void) const {
    return _timeDelay;
} // getTimeDelay


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::sources::PointForce::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains required fields.
    std::string reason = "point source 'PointForce'.";
    size_t numRequired = 0;
    const size_t maxRequired = 2;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "displacement";

    switch (_formulation) {
    case QUASISTATIC:
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        reason = "point source 'PointForce' with dynamic formulation";
        requiredFields[numRequired++] = "velocity";
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::sources::PointForce::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());
    integrator->createLabelDS(solution, solution.getMesh().getDimension());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
std::vector<pylith::feassemble::Constraint*>
pylith::sources::PointForce::createConstraints(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createConstraints(solution="<<solution.getLabel()<<")");

    std::vector<pylith::feassemble::Constraint*> constraints;

    PYLITH_METHOD_RETURN(constraints);
} // createConstraints


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::sources::PointForce::createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("auxiliary field");

    assert(_auxiliaryFactory);
    assert(_scales);
    _auxiliaryFactory->initialize(auxiliaryField, *_scales, domainMesh.getDimension());

    // Add auxiliary subfields
    _auxiliaryFactory->addSourceLocation(_location.data(), _location.size());
    _auxiliaryFactory->addMomentTensor(_momentTensor.data(), _momentTensor.size());
    _auxiliaryFactory->addMagnitude(_magnitude);
    _auxiliaryFactory->addOriginTime(_originTime);
    _auxiliaryFactory->addDominantFrequency(_dominantFrequency);
    _auxiliaryFactory->addTimeDelay(_timeDelay);

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    _auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::sources::PointForce::createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::sources::PointForce::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                  const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    // Nothing to update for now - source time function is evaluated in kernels

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::sources::PointForce::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::sources::PointForce::_updateKernelConstants(const PylithReal dt) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateKernelConstants(dt="<<dt<<")");

    // Store source parameters as kernel constants
    const size_t numConstants = 4;
    _kernelConstants.resize(numConstants);
    _kernelConstants[0] = _originTime;
    _kernelConstants[1] = _dominantFrequency;
    _kernelConstants[2] = _timeDelay;
    _kernelConstants[3] = _magnitude;

    PYLITH_METHOD_END;
} // _updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::sources::PointForce::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                 const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    const int spaceDim = coordsys->getSpaceDim();

    std::vector<ResidualKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        // For quasistatic, point force contributes to displacement equation
        PetscPointFn* f0u = (spaceDim == 2) ?
                           pylith::fekernels::PointForce2D::g0v_pointforce :
                           pylith::fekernels::PointForce3D::g0v_pointforce;
        PetscPointFn* f1u = NULL;

        kernels.resize(1);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        break;
    } // QUASISTATIC
    case DYNAMIC:
    case DYNAMIC_IMEX: {
        // For dynamic, point force contributes to velocity equation (RHS)
        PetscPointFn* g0v = (spaceDim == 2) ?
                           pylith::fekernels::PointForce2D::g0v_pointforce :
                           pylith::fekernels::PointForce3D::g0v_pointforce;
        PetscPointFn* g1v = NULL;

        kernels.resize(1);
        kernels[0] = ResidualKernels("velocity", pylith::feassemble::Integrator::RHS, g0v, g1v);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::sources::PointForce::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                 const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    // Point force does not contribute to Jacobian (no dependence on solution)
    std::vector<JacobianKernels> kernels;

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// End of file
