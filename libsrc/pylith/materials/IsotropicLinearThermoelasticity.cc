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

#include "pylith/materials/IsotropicLinearThermoelasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryThermoelasticity.hh" // USES AuxiliaryFactoryThermoelasticity
#include "pylith/fekernels/IsotropicLinearThermoelasticity.hh" // USES IsotropicLinearThermoelasticity kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearThermoelasticity::IsotropicLinearThermoelasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryThermoelasticity),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearthermoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearThermoelasticity::~IsotropicLinearThermoelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearThermoelasticity::deallocate(void) {
    RheologyThermoelasticity::deallocate();
    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearThermoelasticity::useReferenceState(const bool value) {
    PYLITH_JOURNAL_DEBUG("useReferenceState(value="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearThermoelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryThermoelasticity*
pylith::materials::IsotropicLinearThermoelasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearThermoelasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addAuxiliarySubfields(void)");

    assert(_auxiliaryFactory);

    // Add rheology-specific subfields (shear modulus, bulk modulus)
    // These are queried via vs/vp from the database
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual.
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelf1u_implicit(coordsys="<<coordsys<<")");

    assert(coordsys);
    const int spaceDim = coordsys->getSpaceDim();

    PetscPointFunc kernel = NULL;
    switch (spaceDim) {
    case 2:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::f1u_PlaneStrain;
        break;
    case 3:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::f1u_3D;
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(kernel);
} // getKernelf1u_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual (explicit).
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelg1u_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelg1u_explicit(coordsys="<<coordsys<<")");

    // For explicit, we need the negative of the implicit stress kernel
    // In practice, this is often handled by negating in the integrator
    // For now, return NULL and handle in the material class
    PYLITH_METHOD_RETURN(NULL);
} // getKernelg1u_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian.
PetscPointJac
pylith::materials::IsotropicLinearThermoelasticity::getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelJf3uu(coordsys="<<coordsys<<")");

    assert(coordsys);
    const int spaceDim = coordsys->getSpaceDim();

    PetscPointJac kernel = NULL;
    switch (spaceDim) {
    case 2:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::Jf3uu_PlaneStrain;
        break;
    case 3:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::Jf3uu_3D;
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(kernel);
} // getKernelJf3uu


// ---------------------------------------------------------------------------------------------------------------------
// Get coupling kernel Jf2_uT for stress-temperature coupling.
PetscPointJac
pylith::materials::IsotropicLinearThermoelasticity::getKernelJf2uT(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelJf2uT(coordsys="<<coordsys<<")");

    assert(coordsys);
    const int spaceDim = coordsys->getSpaceDim();

    PetscPointJac kernel = NULL;
    switch (spaceDim) {
    case 2:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::Jf2uT_PlaneStrain;
        break;
    case 3:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::Jf2uT_3D;
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(kernel);
} // getKernelJf2uT


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for LHS residual (implicit).
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelf1T_implicit(coordsys="<<coordsys<<")");

    PYLITH_METHOD_RETURN(pylith::fekernels::IsotropicLinearThermoelasticity::f1T);
} // getKernelf1T_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for RHS residual (explicit).
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelg1T_explicit(coordsys="<<coordsys<<")");

    // Return NULL for now - explicit time integration not fully implemented
    PYLITH_METHOD_RETURN(NULL);
} // getKernelg1T_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Get thermal conductivity kernel for LHS Jacobian.
PetscPointJac
pylith::materials::IsotropicLinearThermoelasticity::getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelJf3TT(coordsys="<<coordsys<<")");

    PYLITH_METHOD_RETURN(pylith::fekernels::IsotropicLinearThermoelasticity::Jf3TT);
} // getKernelJf3TT


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelCauchyStressVector(coordsys="<<coordsys<<")");

    assert(coordsys);
    const int spaceDim = coordsys->getSpaceDim();

    PetscPointFunc kernel = NULL;
    switch (spaceDim) {
    case 2:
        kernel = pylith::fekernels::IsotropicLinearThermoelasticity::cauchyStress_PlaneStrain;
        break;
    case 3:
        // Would need to implement cauchyStress_3D
        kernel = NULL;
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearThermoelasticity::getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getKernelHeatFluxVector(coordsys="<<coordsys<<")");

    PYLITH_METHOD_RETURN(pylith::fekernels::IsotropicLinearThermoelasticity::heatFlux_asVector);
} // getKernelHeatFluxVector


// End of file
