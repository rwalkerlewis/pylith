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

#include "pylith/materials/Heat.hh" // implementation of object methods

#include "pylith/materials/RheologyHeat.hh" // HASA RheologyHeat
#include "pylith/materials/AuxiliaryFactoryHeat.hh" // USES AuxiliaryFactoryHeat
#include "pylith/materials/DerivedFactoryHeat.hh" // USES DerivedFactoryHeat
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/HeatEquation.hh" // USES HeatEquation kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/scales/Scales.hh" // USES Scales

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Heat::Heat(void) :
    _useHeatSource(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryHeat) {
    pylith::utils::PyreComponent::setName("heat");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Heat::~Heat(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Heat::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include heat source?
void
pylith::materials::Heat::useHeatSource(const bool value) {
    PYLITH_COMPONENT_DEBUG("useHeatSource(value="<<value<<")");

    _useHeatSource = value;
} // useHeatSource


// ---------------------------------------------------------------------------------------------------------------------
// Include heat source?
bool
pylith::materials::Heat::useHeatSource(void) const {
    return _useHeatSource;
} // useHeatSource


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Heat::setBulkRheology(pylith::materials::RheologyHeat* const rheology) {
    PYLITH_COMPONENT_DEBUG("setBulkRheology(rheology="<<rheology<<")");

    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyHeat*
pylith::materials::Heat::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Heat::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains required fields.
    std::string reason = "material 'Heat'.";
    size_t numRequired = 0;
    const size_t maxRequired = 2;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "temperature";

    switch (_formulation) {
    case QUASISTATIC:
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        break;
    } // switch
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Heat::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());
    integrator->createLabelDS(solution, solution.getMesh().getDimension());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);
    _setKernelsUpdateStateVars(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::Heat::createAuxiliaryField(const pylith::topology::Field& solution,
                                              const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryHeat* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_scales);
    auxiliaryFactory->initialize(auxiliaryField, *_scales, domainMesh.getDimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // ---------------------------------
    // Required Auxiliary
    auxiliaryFactory->addDensity(); // 0 Density
    auxiliaryFactory->addSpecificHeat(); // 1 Specific Heat
    auxiliaryFactory->addThermalConductivity(); // 2 Thermal Conductivity

    // ---------------------------------
    // Optional Auxiliary
    if (_useHeatSource) {
        auxiliaryFactory->addHeatSource();
    } // if

    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    pythia::journal::debug_t debug("heat.view_auxiliary_field");
    if (debug.state()) {
        auxiliaryField->view("Heat auxiliary field");
    } // if

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::Heat::createDerivedField(const pylith::topology::Field& solution,
                                            const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 1) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->setLabel("derived field");

    assert(_scales);
    _derivedFactory->initialize(derivedField, *_scales, domainMesh.getDimension());
    _derivedFactory->addSubfields();

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->createOutputVector();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Get default PETSc solver options appropriate for material.
pylith::utils::PetscOptions*
pylith::materials::Heat::getSolverDefaults(const bool isParallel,
                                           const bool hasFault) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSolverDefaults(isParallel="<<isParallel<<", hasFault="<<hasFault<<")");

    pylith::utils::PetscOptions* options = new pylith::utils::PetscOptions();assert(options);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        options->add("-ts_type", "beuler");
        options->add("-ts_exact_final_time", "matchstep");

        if (!isParallel) {
            options->add("-pc_type", "lu");
        } else {
            options->add("-pc_type", "ml");
            options->add("-ksp_type", "gmres");
        } // if/else

        break;
    case pylith::problems::Physics::DYNAMIC:
    case pylith::problems::Physics::DYNAMIC_IMEX:
        options->add("-ts_type", "beuler");
        options->add("-ts_exact_final_time", "matchstep");

        if (!isParallel) {
            options->add("-pc_type", "lu");
        } else {
            options->add("-pc_type", "ml");
            options->add("-ksp_type", "gmres");
        } // if/else
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '" << _formulation << "'.");
    } // switch

    PYLITH_METHOD_RETURN(options);
} // getSolverDefaults


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Heat::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Heat::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::Heat::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ----------------------------------------------------------------------
// Set kernels for residual.
void
pylith::materials::Heat::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                             const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        // Temperature
        PetscPointFn* f0T = _useHeatSource ? pylith::fekernels::HeatEquation::f0T_source : NULL;
        PetscPointFn* f1T = _rheology->getKernelf1T_implicit(coordsys);

        kernels.resize(1);
        kernels[0] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        // Temperature (time derivative on LHS, heat flux on RHS)
        PetscPointFn* f0T = pylith::fekernels::HeatEquation::f0T_timedep;
        PetscPointFn* f1T = NULL;
        PetscPointFn* g0T = _useHeatSource ? pylith::fekernels::HeatEquation::g0T_source : NULL;
        PetscPointFn* g1T = _rheology->getKernelg1T_explicit(coordsys);

        kernels.resize(2);
        kernels[0] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
        kernels[1] = ResidualKernels("temperature", pylith::feassemble::Integrator::RHS, g0T, g1T);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    // Add any MMS body force kernels.
    kernels.insert(kernels.end(), _mmsBodyForceKernels.begin(), _mmsBodyForceKernels.end());

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ----------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::Heat::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                             const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<",solution="<<solution.getLabel()<<")");
    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();

    std::vector<JacobianKernels> kernels;

    integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

    switch (_formulation) {
    case QUASISTATIC: {
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS;

        PetscPointJacFn* Jf0TT = NULL;
        PetscPointJacFn* Jf1TT = NULL;
        PetscPointJacFn* Jf2TT = NULL;
        PetscPointJacFn* Jf3TT = _rheology->getKernelJf3TT(coordsys);

        kernels.resize(1);
        kernels[0] = JacobianKernels("temperature", "temperature", equationPart, Jf0TT, Jf1TT, Jf2TT, Jf3TT);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    case DYNAMIC: {
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS;

        PetscPointJacFn* Jf0TT = pylith::fekernels::HeatEquation::Jf0TT_timedep;
        PetscPointJacFn* Jf1TT = NULL;
        PetscPointJacFn* Jf2TT = NULL;
        PetscPointJacFn* Jf3TT = NULL;

        kernels.resize(1);
        kernels[0] = JacobianKernels("temperature", "temperature", equationPart, Jf0TT, Jf1TT, Jf2TT, Jf3TT);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::Heat::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                    const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;

    // No state variables to update for basic heat equation

    integrator->setKernelsUpdateStateVars(kernels);

    PYLITH_METHOD_END;
} // _setKernelsUpdateStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Heat::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                 const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(1);
    kernels[0] = ProjectKernels("heat_flux", _rheology->getKernelHeatFluxVector(coordsys));

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
