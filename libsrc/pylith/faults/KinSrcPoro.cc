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

#include "KinSrcPoro.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDisretization()

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcPoro::KinSrcPoro(void) :
    _auxiliaryFactory(new pylith::faults::KinSrcAuxiliaryFactory),
    _thicknessFnKernel(NULL),
    _porosityFnKernel(NULL),
    _beta_pFnKernel(NULL),
    _beta_sigmaFnKernel(NULL),
    _permeability_tangentialFnKernel(NULL),
    _permeability_normalFnKernel(NULL),
    _fluid_viscosityFnKernel(NULL),
    _bulk_modulus_negativeFnKernel(NULL),
    _shear_modulus_negativeFnKernel(NULL),
    _bulk_modulus_positiveFnKernel(NULL),
    _shear_modulus_positiveFnKernel(NULL),
    _slipFnKernel(NULL),
    _slipRateFnKernel(NULL),
    _slipAccFnKernel(NULL),
    _auxiliaryField(NULL),
    _originTime(0.0) {}


// ----------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcPoro::~KinSrcPoro(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::faults::KinSrcPoro::deallocate(void) {
    delete _auxiliaryField;_auxiliaryField = NULL;
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ----------------------------------------------------------------------
// Set origin time for earthquake source.
void
pylith::faults::KinSrcPoro::originTime(const PylithReal value) {
    _originTime = value;
} // originTime


// ----------------------------------------------------------------------
// Get origin time for earthquake source.
PylithReal
pylith::faults::KinSrcPoro::originTime(void) const {
    return _originTime;
} // originTime


// ----------------------------------------------------------------------
// Get auxiliary field.
const pylith::topology::Field&
pylith::faults::KinSrcPoro::auxField(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_auxiliaryField);

    PYLITH_METHOD_RETURN(*_auxiliaryField);
} // auxField


// ----------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::faults::KinSrcPoro::auxFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("auxFieldDB(value="<<typeid(value).name()<<")");

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setQueryDB(value);

    PYLITH_METHOD_END;
} // auxFieldDB


// ----------------------------------------------------------------------
// Initialize kinematic (prescribed slip) earthquake source.
void
pylith::faults::KinSrcPoro::initialize(const pylith::topology::Field& faultAuxField,
                                       const spatialdata::units::Nondimensional& normalizer,
                                       const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("initialize(faultAuxField"<<faultAuxField.getLabel()<<", normalizer, cs="<<typeid(cs).name()<<")");

    // Set default discretization of auxiliary subfields to match slip/slip_rate subfield in integrator auxiliary field.
    assert(_auxiliaryFactory);
    const char* slipFieldName = faultAuxField.hasSubfield("slip") ? "slip" : "slip_rate";
    const pylith::topology::FieldBase::Discretization& discretization = faultAuxField.getSubfieldInfo(slipFieldName).fe;
    _auxiliaryFactory->setSubfieldDiscretization("default", discretization.basisOrder, discretization.quadOrder,
                                                 discretization.dimension, discretization.isFaultOnly,
                                                 discretization.cellBasis, discretization.feSpace,
                                                 discretization.isBasisContinuous);

    delete _auxiliaryField;_auxiliaryField = new pylith::topology::Field(faultAuxField.getMesh());assert(_auxiliaryField);
    _auxiliaryField->setLabel("KinSrcPoro auxiliary");
    _auxiliaryFieldSetup(normalizer, cs);
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(faultAuxField, *_auxiliaryField);
    _auxiliaryField->allocate();

    _auxiliaryFactory->setValuesFromDB();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying kinematic earthquake source auxiliary field");
        _auxiliaryField->view("KinSrcPoro auxiliary field");
    } // if

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Set slip values at time t.
// Also set rest of auxiliary fields
void
pylith::faults::KinSrcPoro::updateSlip(PetscVec slipLocalVec,
                                       pylith::topology::Field* faultAuxiliaryField,
                                       const PylithScalar t,
                                       const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlip(slipLocalVec="<<slipLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                     <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip vector.

    PetscPointFunc subfieldKernels[12];
    subfieldKernels[0] = _thicknessFnKernel; // 0
    subfieldKernels[1] = _porosityFnKernel; // 1
    subfieldKernels[2] = _beta_pFnKernel; // 2
    subfieldKernels[3] = _beta_sigmaFnKernel; // 3
    subfieldKernels[4] = _permeability_tangentialFnKernel; // 4
    subfieldKernels[5] = _permeability_normalFnKernel; // 5
    subfieldKernels[6] = _fluid_viscosityFnKernel; // 6
    subfieldKernels[7] = _bulk_modulus_negativeFnKernel; // 7
    subfieldKernels[8] = _shear_modulus_negativeFnKernel; // 8
    subfieldKernels[9] = _bulk_modulus_positiveFnKernel; // 9
    subfieldKernels[10] = _shear_modulus_positiveFnKernel; // 10
    subfieldKernels[11] = _slipFnKernel; // 11

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipLocalVec, subfieldKernels, INSERT_VALUES,
                              slipLocalVec);PYLITH_CHECK_ERROR(err);

    // DEBUG LINES
    VecView(slipLocalVec, PETSC_VIEWER_STDOUT_SELF);

    PYLITH_METHOD_END;
} // updateSlip


// ----------------------------------------------------------------------
// Set slip rate values at time t.
void
pylith::faults::KinSrcPoro::updateSlipRate(PetscVec slipRateLocalVec,
                                           pylith::topology::Field* faultAuxiliaryField,
                                           const PylithScalar t,
                                           const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlipRate(slipRateLocalVec="<<slipRateLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                             <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipRateFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipRateLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip rate vector.

    PetscPointFunc subfieldKernels[12];
    subfieldKernels[0] = _thicknessFnKernel; // 0
    subfieldKernels[1] = _porosityFnKernel; // 1
    subfieldKernels[2] = _beta_pFnKernel; // 2
    subfieldKernels[3] = _beta_sigmaFnKernel; // 3
    subfieldKernels[4] = _permeability_tangentialFnKernel; // 4
    subfieldKernels[5] = _permeability_normalFnKernel; // 5
    subfieldKernels[6] = _fluid_viscosityFnKernel; // 6
    subfieldKernels[7] = _bulk_modulus_negativeFnKernel; // 7
    subfieldKernels[8] = _shear_modulus_negativeFnKernel; // 8
    subfieldKernels[9] = _bulk_modulus_positiveFnKernel; // 9
    subfieldKernels[10] = _shear_modulus_positiveFnKernel; // 10
    subfieldKernels[11] = _slipRateFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    PetscDMLabel dmLabel = NULL;
    PetscInt labelValue = 0;
    err = DMSetAuxiliaryVec(faultAuxiliaryDM, dmLabel, labelValue,
                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipRateLocalVec, subfieldKernels, INSERT_VALUES,
                              slipRateLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // updateSlipRate


// ----------------------------------------------------------------------
// Set slip acceleration values at time t.
void
pylith::faults::KinSrcPoro::updateSlipAcc(PetscVec slipAccLocalVec,
                                          pylith::topology::Field* faultAuxiliaryField,
                                          const PylithScalar t,
                                          const PylithScalar timeScale) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateSlipAcc(slipAccLocalVec="<<slipAccLocalVec<<", faultAuxiliaryField="<<faultAuxiliaryField
                                                           <<", t="<<t<<", timeScale="<<timeScale<<")");

    if (!_slipAccFnKernel || (t < _originTime)) {
        PYLITH_METHOD_END;
    } // if

    assert(slipAccLocalVec);
    assert(_auxiliaryField);

    _setFEConstants(*faultAuxiliaryField); // Constants are attached to the auxiliary field for the slip rate vector.

    PetscPointFunc subfieldKernels[12];
    subfieldKernels[0] = _thicknessFnKernel; // 0
    subfieldKernels[1] = _porosityFnKernel; // 1
    subfieldKernels[2] = _beta_pFnKernel; // 2
    subfieldKernels[3] = _beta_sigmaFnKernel; // 3
    subfieldKernels[4] = _permeability_tangentialFnKernel; // 4
    subfieldKernels[5] = _permeability_normalFnKernel; // 5
    subfieldKernels[6] = _fluid_viscosityFnKernel; // 6
    subfieldKernels[7] = _bulk_modulus_negativeFnKernel; // 7
    subfieldKernels[8] = _shear_modulus_negativeFnKernel; // 8
    subfieldKernels[9] = _bulk_modulus_positiveFnKernel; // 9
    subfieldKernels[10] = _shear_modulus_positiveFnKernel; // 10
    subfieldKernels[11] = _slipAccFnKernel;

    // Create local vector for slip for this source.
    PetscErrorCode err = 0;
    PetscDM faultAuxiliaryDM = faultAuxiliaryField->getDM();
    err = PetscObjectCompose((PetscObject) faultAuxiliaryDM, "dmAux",
                             (PetscObject) _auxiliaryField->getDM());PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) faultAuxiliaryDM, "A",
                             (PetscObject) _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(faultAuxiliaryDM, t, slipAccLocalVec, subfieldKernels, INSERT_VALUES,
                              slipAccLocalVec);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // updateSlipAcc


// TO DO:
// Check whether one needs to call this _setFEConstants in every update function above
// Check if numConstants should be one or not
// ----------------------------------------------------------------------
// Set constants used in finite-element integrations.
void
pylith::faults::KinSrcPoro::_setFEConstants(const pylith::topology::Field& faultAuxField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEConstants(faultAuxField="<<faultAuxField.getLabel()<<")");

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.
    PetscDS prob = NULL;
    PetscDM dmAux = faultAuxField.getDM();assert(dmAux);
    PetscErrorCode err = DMGetDS(dmAux, &prob);PYLITH_CHECK_ERROR(err);assert(prob);

    // Pointwise functions have been set in DS
    const int numConstants = 1;
    PylithScalar constants[numConstants];
    constants[0] = _originTime;
    err = PetscDSSetConstants(prob, numConstants, constants);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _setFEConstants


// End of file
