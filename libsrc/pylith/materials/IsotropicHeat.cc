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

#include "pylith/materials/IsotropicHeat.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryHeat.hh" // USES AuxiliaryFactoryHeat
#include "pylith/fekernels/IsotropicHeat.hh" // USES IsotropicHeat kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicHeat::IsotropicHeat(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryHeat) {
    pylith::utils::PyreComponent::setName("isotropicheat");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicHeat::~IsotropicHeat(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicHeat::deallocate(void) {
    RheologyHeat::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryHeat*
pylith::materials::IsotropicHeat::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicHeat::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // No additional rheology-specific subfields for isotropic heat.
    // The basic fields (density, specific_heat, thermal_conductivity) are added in Heat::createAuxiliaryField

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for LHS residual, F(t,s,\dot{s})
PetscPointFn*
pylith::materials::IsotropicHeat::getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1T_implicit(coordsys="<<coordsys<<")");

    PetscPointFn* f1T = pylith::fekernels::IsotropicHeat::f1T;

    PYLITH_METHOD_RETURN(f1T);
} // getKernelf1T_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for RHS residual, G(t,s).
PetscPointFn*
pylith::materials::IsotropicHeat::getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1T_explicit(coordsys="<<coordsys<<")");

    PetscPointFn* g1T = pylith::fekernels::IsotropicHeat::g1T;

    PYLITH_METHOD_RETURN(g1T);
} // getKernelg1T_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Get thermal conductivity kernel for LHS Jacobian.
PetscPointJacFn*
pylith::materials::IsotropicHeat::getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3TT(coordsys="<<coordsys<<")");

    PetscPointJacFn* Jf3TT = pylith::fekernels::IsotropicHeat::Jf3TT;

    PYLITH_METHOD_RETURN(Jf3TT);
} // getKernelJf3TT


// ---------------------------------------------------------------------------------------------------------------------
// Get heat flux kernel for derived field.
PetscPointFn*
pylith::materials::IsotropicHeat::getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelHeatFluxVector(coordsys="<<coordsys<<")");

    PetscPointFn* heatFlux = pylith::fekernels::IsotropicHeat::heatFlux_asVector;

    PYLITH_METHOD_RETURN(heatFlux);
} // getKernelHeatFluxVector


// End of file
