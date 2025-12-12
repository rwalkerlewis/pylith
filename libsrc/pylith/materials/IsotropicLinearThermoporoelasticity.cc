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

#include "pylith/materials/IsotropicLinearThermoporoelasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryThermoporoelasticity.hh" // USES AuxiliaryFactoryThermoporoelasticity
#include "pylith/fekernels/IsotropicLinearThermoporoelasticity.hh" // USES IsotropicLinearThermoporoelasticity kernels
#include "pylith/fekernels/Thermoporoelasticity.hh" // USES Thermoporoelasticity kernels

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::fekernels::IsotropicLinearThermoporoelasticity IsotropicLinearThermoporoelasticityKernels;
typedef pylith::fekernels::Thermoporoelasticity ThermoporoelasticityKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearThermoporoelasticity::IsotropicLinearThermoporoelasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryThermoporoelasticity),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearthermoporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearThermoporoelasticity::~IsotropicLinearThermoporoelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearThermoporoelasticity::deallocate(void) {
    RheologyThermoporoelasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearThermoporoelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState(value="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearThermoporoelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryThermoporoelasticity*
pylith::materials::IsotropicLinearThermoporoelasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearThermoporoelasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // Poroelastic fields
    _auxiliaryFactory->addBiotCoefficient();
    _auxiliaryFactory->addBiotModulus();
    _auxiliaryFactory->addDrainedBulkModulus();
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addIsotropicPermeability();

    // Thermal fields
    _auxiliaryFactory->addReferenceTemperature();
    _auxiliaryFactory->addThermalExpansionCoeff();
    _auxiliaryFactory->addFluidThermalExpansion();
    _auxiliaryFactory->addThermalConductivity();
    _auxiliaryFactory->addSpecificHeat();

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ============================= RHS (Explicit) ==================================== //

// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelg1u_explicit(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1u_explicit(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* g1u = NULL;

    switch (spaceDim) {
    case 2:
        g1u = IsotropicLinearThermoporoelasticityKernels::f1u_PlaneStrain;
        break;
    case 3:
        g1u = IsotropicLinearThermoporoelasticityKernels::f1u_3D;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(g1u);
} // getKernelg1u_explicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelg1p_explicit(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1p_explicit(coordsys="<<typeid(coordsys).name()<<", gravityField="<<gravityField<<")");

    PetscPointFn* g1p = NULL;

    if (gravityField) {
        g1p = IsotropicLinearThermoporoelasticityKernels::f1p_darcy_grav;
    } else {
        g1p = IsotropicLinearThermoporoelasticityKernels::f1p_darcy;
    }

    PYLITH_METHOD_RETURN(g1p);
} // getKernelg1p_explicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelg1T_explicit(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1T_explicit(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::f1T_heatflux);
} // getKernelg1T_explicit


// =============================== LHS (Implicit) =================================== //

// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelf1u_implicit(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1u_implicit(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* f1u = NULL;

    switch (spaceDim) {
    case 2:
        f1u = IsotropicLinearThermoporoelasticityKernels::f1u_PlaneStrain;
        break;
    case 3:
        f1u = IsotropicLinearThermoporoelasticityKernels::f1u_3D;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1u_implicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelf0p_implicit(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p_implicit(coordsys="<<typeid(coordsys).name()<<", useSourceDensity="<<useSourceDensity<<")");

    PetscPointFn* f0p = NULL;

    if (useSourceDensity) {
        f0p = IsotropicLinearThermoporoelasticityKernels::f0p_source;
    } else {
        f0p = IsotropicLinearThermoporoelasticityKernels::f0p_fluidContent;
    }

    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_implicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelf1p_implicit(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1p_implicit(coordsys="<<typeid(coordsys).name()<<", gravityField="<<gravityField<<")");

    PetscPointFn* f1p = NULL;

    if (gravityField) {
        f1p = IsotropicLinearThermoporoelasticityKernels::f1p_darcy_grav;
    } else {
        f1p = IsotropicLinearThermoporoelasticityKernels::f1p_darcy;
    }

    PYLITH_METHOD_RETURN(f1p);
} // getKernelf1p_implicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelf0T_implicit(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool useHeatSource) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0T_implicit(coordsys="<<typeid(coordsys).name()<<", useHeatSource="<<useHeatSource<<")");

    PetscPointFn* f0T = NULL;

    if (useHeatSource) {
        f0T = ThermoporoelasticityKernels::f0T_source;
    } else {
        f0T = ThermoporoelasticityKernels::f0T_transient;
    }

    PYLITH_METHOD_RETURN(f0T);
} // getKernelf0T_implicit


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelf1T_implicit(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1T_implicit(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::f1T_heatflux);
} // getKernelf1T_implicit


// ============================= LHS Jacobians ==================================== //

// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf3uu(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3uu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3uu = NULL;

    switch (spaceDim) {
    case 2:
        Jf3uu = IsotropicLinearThermoporoelasticityKernels::Jf3uu_PlaneStrain;
        break;
    case 3:
        Jf3uu = IsotropicLinearThermoporoelasticityKernels::Jf3uu_3D;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJf3uu


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf2up(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2up(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf2up);
} // getKernelJf2up


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf2uT(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2uT(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf2uT);
} // getKernelJf2uT


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf0pp(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pp(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf0pp);
} // getKernelJf0pp


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf0pe(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pe(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf0pe);
} // getKernelJf0pe


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf0pT(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pT(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf0pT);
} // getKernelJf0pT


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf3pp(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3pp(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf3pp);
} // getKernelJf3pp


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf0TT(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0TT(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf0TT);
} // getKernelJf0TT


// ---------------------------------------------------------------------------------------------------------------------
PetscPointJacFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelJf3TT(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3TT(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::Jf3TT);
} // getKernelJf3TT


// ============================ DERIVED FIELDS ========================== //

// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelCauchyStressVector(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel = NULL;

    switch (spaceDim) {
    case 2:
        kernel = IsotropicLinearThermoporoelasticityKernels::cauchyStress_PlaneStrain;
        break;
    case 3:
        // TODO: Add 3D version
        kernel = IsotropicLinearThermoporoelasticityKernels::cauchyStress_PlaneStrain;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown spatial dimension ("<<spaceDim<<").");
    } // switch

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelFluidContent(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelFluidContent(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::fluidContent);
} // getKernelFluidContent


// ---------------------------------------------------------------------------------------------------------------------
PetscPointFn*
pylith::materials::IsotropicLinearThermoporoelasticity::getKernelHeatFluxVector(
    const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelHeatFluxVector(coordsys="<<typeid(coordsys).name()<<")");

    PYLITH_METHOD_RETURN(IsotropicLinearThermoporoelasticityKernels::heatFlux);
} // getKernelHeatFluxVector


// End of file
