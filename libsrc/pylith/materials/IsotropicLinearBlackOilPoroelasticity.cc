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

#include "pylith/materials/IsotropicLinearBlackOilPoroelasticity.hh"
#include "pylith/materials/AuxiliaryFactoryBlackOilPoroelastic.hh"

#include "pylith/fekernels/IsotropicLinearBlackOilPoroelasticity.hh"
#include "pylith/fekernels/Elasticity.hh"

#include "pylith/utils/journals.hh"
#include "pylith/utils/error.hh"

#include "spatialdata/geocoords/CoordSys.hh"

#include <typeinfo>

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearBlackOilPoroelasticity::IsotropicLinearBlackOilPoroelasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryBlackOilPoroelastic),
    _useReferenceState(false),
    _useTensorPermeability(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearblackoilporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearBlackOilPoroelasticity::~IsotropicLinearBlackOilPoroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::deallocate(void) {
    RheologyPoroelasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearBlackOilPoroelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::useTensorPermeability(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTensorPermeability="<<value<<")");

    _useTensorPermeability = value;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
bool
pylith::materials::IsotropicLinearBlackOilPoroelasticity::useTensorPermeability(void) const {
    return _useTensorPermeability;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryPoroelastic*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addDrainedBulkModulus();
    _auxiliaryFactory->addBiotCoefficient();
    _auxiliaryFactory->addBiotModulus();

    // Black oil specific fields
    _auxiliaryFactory->addReferencePressure();
    _auxiliaryFactory->addFluidCompressibility();
    _auxiliaryFactory->addFluidCompressibilityCoefficient();
    _auxiliaryFactory->addViscosityCoefficient();

    if (_useTensorPermeability) {
        _auxiliaryFactory->addTensorPermeability();
    } else {
        _auxiliaryFactory->addIsotropicPermeability();
    }

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ================================ RHS ========================================

// ---------------------------------------------------------------------------------------------------------------------
// Select g0p function.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                                                                       const bool _useBodyForce,
                                                                       const bool _gravityField,
                                                                       const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg0p="<<typeid(coordsys).name()<<")");

    // For now, return NULL as this is only used for the dynamic case
    // TODO: Implement dynamic kernels for black oil model
    PetscPointFn* g0p = NULL;

    PYLITH_METHOD_RETURN(g0p);
} // getKernelg0p


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy velocity kernel for explicit time stepping
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                                const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1p_explicit="<<typeid(coordsys).name()<<")");

    // For now, return NULL as this is only used for explicit time stepping
    // TODO: Implement explicit kernels for black oil model
    PetscPointFn* g1p = NULL;

    PYLITH_METHOD_RETURN(g1p);
} // getKernelg1p_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelg1v_explicit(coordsys="<<typeid(coordsys).name()<<")");

    // For now, return NULL as this is only used for explicit time stepping
    PetscPointFn* g1v = NULL;

    PYLITH_METHOD_RETURN(g1v);
} // getKernelg1v_explicit


// =============================== LHS =========================================

// ---------------------------------------------------------------------------------------------------------------------
// Get variation in fluid content kernel for LHS residual, explicit time stepping
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p_explicit="<<typeid(coordsys).name()<<")");

    // For now, return NULL as this is only used for explicit time stepping
    PetscPointFn* f0p = NULL;

    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_explicit


// ---------------------------------------------------------------------------------------------------------------------
// Select implicit f0p function.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                                const bool _useBodyForce,
                                                                                const bool _gravityField,
                                                                                const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitSourceDensity = _useSourceDensity ? 0x1 : 0x0;
    const int bitUse = bitSourceDensity;

    PetscPointFn* f0p = NULL;

    switch (bitUse) {
    case 0x0:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x1:
        f0p = (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::f0p_implicit_source :
              (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::f0p_implicit_source :
              NULL;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ").");
    } // switch

    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1u(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();

    PetscPointFn* f1u =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::f1u :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::f1u :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1u_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy velocity kernel for implicit time stepping
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                                const bool _useBodyForce,
                                                                                const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1p_implicit="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();

    // For black oil, we currently only support isotropic permeability without gravity
    // TODO: Add support for tensor permeability and gravity
    PetscPointFn* f1p =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::f1p :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::f1p :
        NULL;

    PYLITH_METHOD_RETURN(f1p);
} // getKernelf1p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get poroelastic constants kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3uu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf3uu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJf3uu


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2up(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf2up =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf2up :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf2up :
        NULL;

    PYLITH_METHOD_RETURN(Jf2up);
} // getKernelJf2up


// ---------------------------------------------------------------------------------------------------------------------
// Get lambda kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2ue(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf2ue =
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf2ue :
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf2ue :
        NULL;

    PYLITH_METHOD_RETURN(Jf2ue);
} // getKernelJf2ue


// ---------------------------------------------------------------------------------------------------------------------
// Get Specific storage kernel for LHS Jacobian.
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf0pp =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf0pp :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf0pp :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pp);
} // getKernelJf0pp


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy Conductivity kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3pp =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf3pp :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf3pp :
        NULL;

    PYLITH_METHOD_RETURN(Jf3pp);
} // getKernelJf3pp


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian.
PetscPointJacFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pe(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf0pe =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::Jf0pe :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::Jf0pe :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pe);
} // getKernelJf0pe


// =========================== DERIVED FIELDS ==================================

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();

    PetscPointFn* kernel =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::cauchyStress_infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ---------------------------------------------------------------------------------------------------------------------
// Get water content kernel for derived field.
PetscPointFn*
pylith::materials::IsotropicLinearBlackOilPoroelasticity::getKernelWaterContent(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelWaterContent(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel =
        (3 == spaceDim) ?  pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D::waterContent_asScalar :
        (2 == spaceDim) ?  pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain::waterContent_asScalar :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelWaterContent


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                                const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, implicit.
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::addKernelsUpdateStateVarsImplicit(std::vector<ProjectKernels>* kernels,
                                                                                           const spatialdata::geocoords::CoordSys* coordsys,
                                                                                           const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsImplicit(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    // Black oil model does not currently have state variable updates
    // TODO: Add porosity update if needed

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsImplicit


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, explicit.
void
pylith::materials::IsotropicLinearBlackOilPoroelasticity::addKernelsUpdateStateVarsExplicit(std::vector<ProjectKernels>* kernels,
                                                                                           const spatialdata::geocoords::CoordSys* coordsys,
                                                                                           const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsExplicit(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    // Black oil model does not currently have state variable updates

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsExplicit


// End of file
