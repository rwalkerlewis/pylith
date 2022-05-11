// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesivePoroKin.hh" // implementation of object methods

#include "pylith/faults/KinSrc.hh" // USES KinSrc
#include "pylith/faults/AuxiliaryFactoryKinematic.hh" // USES AuxiliaryFactoryKinematic
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/feassemble/ConstraintSimple.hh" // USES ConstraintSimple

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDiscretization()
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/fekernels/FaultCohesivePoroKin.hh" // USES FaultCohesivePoroKin

#include "pylith/utils/EventLogger.hh" // USES EventLogger
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

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorInterface::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorInterface::JacobianKernels JacobianKernels;

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class _FaultCohesivePoroKin {
            // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////
public:

            static const char* pyreComponent;

        };
        const char* _FaultCohesivePoroKin::pyreComponent = "faultcohesiveporokin";

        // _FaultCohesivePoroKin

    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesivePoroKin::FaultCohesivePoroKin(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryKinematicPoro),
    _slipVecRupture(NULL),
    _slipVecTotal(NULL) {
    pylith::utils::PyreComponent::setName(_FaultCohesivePoroKin::pyreComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesivePoroKin::~FaultCohesivePoroKin(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::FaultCohesivePoroKin::deallocate(void) {
    FaultCohesive::deallocate();

    PetscErrorCode err = VecDestroy(&_slipVecRupture);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_slipVecTotal);PYLITH_CHECK_ERROR(err);
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
    _ruptures.clear(); // :TODO: Use shared pointers for earthquake ruptures
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set kinematic earthquake ruptures.
void
pylith::faults::FaultCohesivePoroKin::setEqRuptures(const char* const * names,
                                                const int numNames,
                                                KinSrc** ruptures,
                                                const int numRuptures) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setEqRuptures(names="<<names<<", numNames="<<numNames<<", ruptures="<<ruptures<<", numRuptures="<<numRuptures<<")");

    assert(numNames == numRuptures);

    // :TODO: Use shared pointers for earthquake ruptures
    _ruptures.clear();
    for (int i = 0; i < numRuptures; ++i) {
        if (!ruptures[i]) {
            std::ostringstream msg;
            msg << "Null earthquake rupture object for earthquake rupture '" << names[i] << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _ruptures[std::string(names[i])] = ruptures[i];
    } // for

    PYLITH_METHOD_END;
} // setEqRuptures


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::faults::FaultCohesivePoroKin::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution=" << solution.getLabel() << ")");

    if (!solution.hasSubfield("lagrange_multiplier_fault")) {
        std::ostringstream msg;
        msg << "Cannot find 'lagrange_multiplier_fault' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    /** TO DO **
     * Implement or verify field::hasSubfield("pressure", "trace_strain", "fault_pressure") is implemented
     * Seems that this is fine, at least no need to change the function field::hasSubfield
     */
    // ** TO DO **
    // Seems that pressure is not part of the fault, but part of poroelastic material implementation

    if (!solution.hasSubfield("pressure")) {
        std::ostringstream msg;
        msg << "Cannot find 'pressure' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    if (!solution.hasSubfield("fault_pressure")) {
        std::ostringstream msg;
        msg << "Cannot find 'fault_pressure' subfield in solution field for poroelastic diffusive fault implementation in component '"
            << PyreComponent::getIdentifier() << "'.";
        throw std::runtime_error(msg.str());
    } // if

    switch (_formulation) {
    case QUASISTATIC:
        if (!solution.hasSubfield("displacement")) {
            std::ostringstream msg;
            msg << "Cannot find 'displacement' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
        if (!solution.hasSubfield("trace_strain")) {
            std::ostringstream msg;
            msg << "Cannot find 'trace_strain' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC_IMEX:
        if (!solution.hasSubfield("velocity")) {
            std::ostringstream msg;
            msg << "Cannot find 'velocity' subfield in solution field for fault implementation in component '"
                << PyreComponent::getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::faults::FaultCohesivePoroKin::createAuxiliaryField(const pylith::topology::Field& solution,
                                                       const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution=" << solution.getLabel() << ", domainMesh=)" << typeid(domainMesh).name() << ")");

    assert(_normalizer);

    pylith::topology::Field *auxiliaryField = new pylith::topology::Field(domainMesh);
    assert(auxiliaryField);
    auxiliaryField->setLabel("FaultCohesivePoroKin auxiliary field");

    // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    assert(_auxiliaryFactory);
    const pylith::topology::FieldBase::Discretization &discretization = solution.getSubfieldInfo("lagrange_multiplier_fault").fe;
    bool isFaultOnly = false;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder, -1, isFaultOnly,
                                                 discretization.cellBasis, discretization.feSpace, discretization.isBasisContinuous);

    // // Set default discretization of auxiliary subfields to match lagrange_multiplier_fault subfield in solution.
    // assert(_auxiliaryFactory);
    // const pylith::topology::FieldBase::Discretization &discretization1 =
    // solution.getSubfieldInfo("fault_pressure").fe;
    // isFaultOnly = false;
    // _auxiliaryFactory->setSubfieldDiscretization("default", discretization1.basisOrder, discretization1.quadOrder,
    // -1, isFaultOnly,
    //                                              discretization1.cellBasis, discretization1.feSpace,
    // discretization1.isBasisContinuous);

    // ** TO DO **
    // Set discretization of aux fields to match also
    // pressure (\Gamma^+ and - sides), trace_strain (\Gamma^+ and - sides), fault_pressure (only \Gamma^f)

    assert(_auxiliaryFactory);
    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, solution.getSpaceDim());

    // ** TO DO **
    // implement the following add* functions
    /** -- numA : number of auxiliary fields
     ***** Required fields
     * - 0: thickness(1)
     * - 1: porosity(1)
     * - 2: beta_p(1)
     * - 3: beta_sigma(1)
     * - 4: fault_permeability(dim*2)
     * - 5: fluid_viscosity(1)
     * - numA - 3: body_force(dim)
     * - numA - 2: source(1)
     * - numA - 1: slip(dim)
     */
    _auxiliaryFactory->addThickness(); // 0
    _auxiliaryFactory->addPorosity(); // 1
    _auxiliaryFactory->addBetaP(); // 2
    _auxiliaryFactory->addBetaSigma(); // 3
    _auxiliaryFactory->addFaultPermeability(); // 4
    _auxiliaryFactory->addFluidViscosity(); // 5

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    switch (_formulation) {
    case QUASISTATIC:
        _auxiliaryFactory->addSlip(); // numA - 1
        break;
    case DYNAMIC_IMEX:
        _auxiliaryFactory->addSlipAcceleration(); // numA - 1
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' time-stepping formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    // We don't populate the auxiliary field via a spatial database, because they will be set from the earthquake
    // rupture.

    // Initialize auxiliary fields for kinematic ruptures.
    assert(auxiliaryField);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        KinSrcPoro *src = r_iter->second;
        assert(src);
        src->initialize(*auxiliaryField, *_normalizer, solution.getMesh().getCoordSys());
    } // for

    // Create local PETSc vector to hold current slip.
    PetscErrorCode err = 0;
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_slipVecRupture);
    PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(auxiliaryField->getDM(), &_slipVecTotal);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Update auxiliary fields at beginning of time step.
void
pylith::faults::FaultCohesivePoroKin::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                       const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateAuxiliaryField(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    switch (_formulation) {
    case QUASISTATIC:
        this->_updateSlip(auxiliaryField, t);
        break;
    case DYNAMIC_IMEX:
        this->_updateSlipAcceleration(auxiliaryField, t);
        break;
    case DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::faults::FaultCohesivePoroKin::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // _getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Update slip subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesivePoroKin::_updateSlip(pylith::topology::Field* auxiliaryField,
                                              const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        err = VecSet(_slipVecRupture, 0.0);PYLITH_CHECK_ERROR(err);

        KinSrcPoro* src = r_iter->second;assert(src);
        src->updateSlip(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip");
    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar* slipArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipArray);PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlip = 0; p < pEnd; ++p) {
        const PetscInt slipDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipDof; ++iDof, ++iSlip) {
            auxiliaryArray[slipOff+iDof] = slipArray[iSlip];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipArray);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting slip.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlip


// ------------------------------------------------------------------------------------------------
// Update slip rate subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesivePoroKin::_updateSlipRate(pylith::topology::Field* auxiliaryField,
                                                  const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlipRate(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip rate subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        err = VecSet(_slipVecRupture, 0.0);PYLITH_CHECK_ERROR(err);

        KinSrcPoro* src = r_iter->second;assert(src);
        src->updateSlipRate(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip_rate");
    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar* slipRateArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipRateArray);PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlipRate = 0; p < pEnd; ++p) {
        const PetscInt slipRateDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipRateOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipRateDof; ++iDof, ++iSlipRate) {
            auxiliaryArray[slipRateOff+iDof] = slipRateArray[iSlipRate];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipRateArray);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting slip rate.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlipRate


// ------------------------------------------------------------------------------------------------
// Update slip acceleration subfield in auxiliary field at beginning of time step.
void
pylith::faults::FaultCohesivePoroKin::_updateSlipAcceleration(pylith::topology::Field* auxiliaryField,
                                                          const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_updateSlipAcceleration(auxiliaryField="<<auxiliaryField<<", t="<<t<<")");

    assert(auxiliaryField);
    assert(_normalizer);

    // Update slip acceleration subfield at current time step
    PetscErrorCode err = VecSet(_slipVecTotal, 0.0);PYLITH_CHECK_ERROR(err);
    const srcs_type::const_iterator rupturesEnd = _ruptures.end();
    for (srcs_type::iterator r_iter = _ruptures.begin(); r_iter != rupturesEnd; ++r_iter) {
        err = VecSet(_slipVecRupture, 0.0);PYLITH_CHECK_ERROR(err);

        KinSrcPoro* src = r_iter->second;assert(src);
        src->updateSlipAcc(_slipVecRupture, auxiliaryField, t, _normalizer->getTimeScale());
        err = VecAYPX(_slipVecTotal, 1.0, _slipVecRupture);
    } // for

    // Transfer slip values from local PETSc slip vector to fault auxiliary field.
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField, "slip_acceleration");
    PylithScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    const PylithScalar* slipAccArray = NULL;
    err = VecGetArrayRead(_slipVecTotal, &slipAccArray);PYLITH_CHECK_ERROR(err);

    for (PetscInt p = pStart, iSlipAcc = 0; p < pEnd; ++p) {
        const PetscInt slipAccDof = auxiliaryVisitor.sectionDof(p);
        const PetscInt slipAccOff = auxiliaryVisitor.sectionOffset(p);
        for (PetscInt iDof = 0; iDof < slipAccDof; ++iDof, ++iSlipAcc) {
            auxiliaryArray[slipAccOff+iDof] = slipAccArray[iSlipAcc];
        } // for
    } // for
    err = VecRestoreArrayRead(_slipVecTotal, &slipAccArray);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(pylith::utils::PyreComponent::getName());
    if (debug.state()) {
        auxiliaryField->view("Fault auxiliary field after setting slip acceleration.");
    } // if

    PYLITH_METHOD_END;
} // _updateSlipAcceleration


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::faults::FaultCohesivePoroKin::_setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                                      const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    std::vector<ResidualKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        // Current solution field is
        // [displacement, pressure, trace_strain, lagrange, fault_pressure]

        // Elasticity equation (displacement/velocity).
        const PetscBdPointFunc f0u_neg = pylith::fekernels::FaultCohesivePoroKin::f0u_neg;
        const PetscBdPointFunc f1u_neg = NULL;

        // Elasticity equation (displacement/velocity).
        const PetscBdPointFunc f0u_pos = pylith::fekernels::FaultCohesivePoroKin::f0u_pos;
        const PetscBdPointFunc f1u_pos = NULL;

        // Trace_strain constraint
        // const PetscBdPointFunc f0e = NULL;
        // const PetscBdPointFunc f1e = NULL;

        // Pressure constraint and Fault_pressure constraint
        PetscBdPointFunc f0p_neg = NULL;
        PetscBdPointFunc f0p_pos = NULL;
        PetscBdPointFunc f0p_fault = NULL;
        PetscBdPointFunc f1p_neg = NULL;
        PetscBdPointFunc f1p_pos = NULL;
        PetscBdPointFunc f1p_fault = NULL;

        const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
        const int bitSource = _useSource ? 0x2 : 0x0;
        const int bitUse = bitBodyForce | bitSource;

        switch (bitUse) {
        case 0x0:
            f0p_neg = pylith::fekernels::FaultCohesivePoroKin::f0p_neg;
            f0p_pos = pylith::fekernels::FaultCohesivePoroKin::f0p_pos;
            f0p_fault = pylith::fekernels::FaultCohesivePoroKin::f0p_fault;
            break;
        case 0x1:
            f0p_neg = pylith::fekernels::FaultCohesivePoroKin::f0p_body_neg;
            f0p_pos = pylith::fekernels::FaultCohesivePoroKin::f0p_body_pos;
            f0p_fault = pylith::fekernels::FaultCohesivePoroKin::f0p_fault_body;
            break;
        case 0x2:
            f0p_neg = pylith::fekernels::FaultCohesivePoroKin::f0p_neg;
            f0p_pos = pylith::fekernels::FaultCohesivePoroKin::f0p_pos;
            f0p_fault = pylith::fekernels::FaultCohesivePoroKin::f0p_fault_source;
            break;
        case 0x3:
            f0p_neg = pylith::fekernels::FaultCohesivePoroKin::f0p_body_neg;
            f0p_pos = pylith::fekernels::FaultCohesivePoroKin::f0p_body_pos;
            f0p_fault = pylith::fekernels::FaultCohesivePoroKin::f0p_fault_body_source;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for f0p f0p_fault residual kernels.");
        } // switch

        // Fault slip constraint equation.
        const PetscBdPointFunc f0l = pylith::fekernels::FaultCohesivePoroKin::f0l_u;
        const PetscBdPointFunc f1l = NULL;

        kernels.resize(6);
        kernels[0] = ResidualKernels("displacement", integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                                     f0u_neg, f1u_neg);
        kernels[1] = ResidualKernels("displacement", integrator_t::LHS, integrator_t::POSITIVE_FACE,
                                     f0u_pos, f1u_pos);
        kernels[2] = ResidualKernels("pressure", integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                                     f0p_neg, f1p_neg);
        kernels[3] = ResidualKernels("pressure", integrator_t::LHS, integrator_t::POSITIVE_FACE,
                                     f0p_pos, f1p_pos);
        kernels[4] = ResidualKernels("lagrange_multiplier_fault", integrator_t::LHS, integrator_t::FAULT_FACE,
                                     f0l, f1l);
        kernels[5] = ResidualKernels("fault_pressure", integrator_t::LHS, integrator_t::FAULT_FACE,
                                     f0p_fault, f1p_fault);
        break;
    } // QUASISTATIC

    case pylith::problems::Physics::DYNAMIC_IMEX:
        break;
    // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual

// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::faults::FaultCohesivePoroKin::_setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                                      const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");
    typedef pylith::feassemble::IntegratorInterface integrator_t;

    // Ten non-trivial kernels,
    // Jful_neg, Jful_pos, Jfp_fp_f, Jfp_fl, Jfp_fp;
    // Jfpp_neg, Jfpp_pos, Jfpp_f_neg, Jfpp_f_pos, Jflu;
    std::vector<JacobianKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        // Displacement
        const PetscBdPointJac Jf0ul_neg = pylith::fekernels::FaultCohesivePoroKin::Jf0ul_neg;
        const PetscBdPointJac Jf1ul_neg = NULL;
        const PetscBdPointJac Jf2ul_neg = NULL;
        const PetscBdPointJac Jf3ul_neg = NULL;

        const PetscBdPointJac Jf0ul_pos = pylith::fekernels::FaultCohesivePoroKin::Jf0ul_pos;
        const PetscBdPointJac Jf1ul_pos = NULL;
        const PetscBdPointJac Jf2ul_pos = NULL;
        const PetscBdPointJac Jf3ul_pos = NULL;

        // Pressure

        // Trace Strain

        // Lagrange Multipliers
        const PetscBdPointJac Jf0lu = pylith::fekernels::FaultCohesivePoroKin::Jf0lu;
        const PetscBdPointJac Jf1lu = NULL;
        const PetscBdPointJac Jf2lu = NULL;
        const PetscBdPointJac Jf3lu = NULL;

        // Fault Pressure
        const PetscBdPointJac Jf0p_fp = pylith::fekernels::FaultCohesivePoroKin::Jf0p_fp;
        const PetscBdPointJac Jf1p_fp = NULL;
        const PetscBdPointJac Jf2p_fp = NULL;
        const PetscBdPointJac Jf3p_fp = pylith::fekernels::FaultCohesivePoroKin::Jf3p_fp;

        const PetscBdPointJac Jf0p_fl = pylith::fekernels::FaultCohesivePoroKin::Jf0p_fl;
        const PetscBdPointJac Jf1p_fl = NULL;
        const PetscBdPointJac Jf2p_fl = NULL;
        const PetscBdPointJac Jf3p_fl = NULL;

        const PetscBdPointJac Jf0p_fp_f = pylith::fekernels::FaultCohesivePoroKin::Jf0p_fp_f;
        const PetscBdPointJac Jf1p_fp_f = NULL;
        const PetscBdPointJac Jf2p_fp_f = NULL;
        const PetscBdPointJac Jf3p_fp_f = pylith::fekernels::FaultCohesivePoroKin::Jf3p_fp_f;

        kernels.resize(6);
        const char *nameDisplacement = "displacement";
        const char *nameLagrangeMultiplier = "lagrange_multiplier_fault";
        const char *namePressure = "pressure";
        const char *nameTraceStrain = "trace_strain";
        const char *nameFaultPressure = "fault_pressure";

        kernels[0] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::NEGATIVE_FACE, Jf0ul_neg, Jf1ul_neg, Jf2ul_neg, Jf3ul_neg);
        kernels[1] = JacobianKernels(nameDisplacement, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::POSITIVE_FACE, Jf0ul_pos, Jf1ul_pos, Jf2ul_pos, Jf3ul_pos);
        kernels[2] = JacobianKernels(nameLagrangeMultiplier, nameDisplacement, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0lu, Jf1lu, Jf2lu, Jf3lu);
        kernels[3] = JacobianKernels(nameFaultPressure, namePressure, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0p_fp, Jf1p_fp, Jf2p_fp, Jf3p_fp);
        kernels[4] = JacobianKernels(nameFaultPressure, nameLagrangeMultiplier, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0p_fl, Jf1p_fl, Jf2p_fl, Jf3p_fl);
        kernels[5] = JacobianKernels(nameFaultPressure, nameFaultPressure, integrator_t::LHS,
                                     integrator_t::FAULT_FACE, Jf0p_fp_f, Jf1p_fp_f, Jf2p_fp_f, Jf3p_fp_f);

        break;
    } // QUASISTATIC
    case pylith::problems::Physics::DYNAMIC_IMEX: {
        break;
    } // DYNAMIC_IMEX
    case pylith::problems::Physics::DYNAMIC:
        PYLITH_COMPONENT_LOGICERROR("Fault implementation is incompatible with 'dynamic' formulation. Use 'dynamic_imex'.");
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<"'.");
    } // switch

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// End of file
