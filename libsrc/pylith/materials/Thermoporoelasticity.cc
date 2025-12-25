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

#include "pylith/materials/Thermoporoelasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyThermoporoelasticity.hh" // USES RheologyThermoporoelasticity
#include "pylith/materials/AuxiliaryFactoryThermoporoelasticity.hh" // USES AuxiliaryFactoryThermoporoelasticity
#include "pylith/materials/DerivedFactoryPoroelasticity.hh" // USES DerivedFactoryPoroelasticity
#include "pylith/fekernels/Thermoporoelasticity.hh" // USES Thermoporoelasticity kernels
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity::bulkDensity_asScalar
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity strain kernels
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;
typedef pylith::fekernels::Thermoporoelasticity ThermoporoelasticityKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Thermoporoelasticity::Thermoporoelasticity(void) :
    _useBodyForce(false),
    _useSourceDensity(false),
    _useHeatSource(false),
    _useReferenceState(false),
    _useStateVars(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryPoroelasticity) {
    pylith::utils::PyreComponent::setName("thermoporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Thermoporoelasticity::~Thermoporoelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Thermoporoelasticity::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // Held by Python; Python must deallocate.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::Thermoporoelasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::Thermoporoelasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include source density?
void
pylith::materials::Thermoporoelasticity::useSourceDensity(const bool value) {
    PYLITH_COMPONENT_DEBUG("useSourceDensity(value="<<value<<")");

    _useSourceDensity = value;
} // useSourceDensity


// ---------------------------------------------------------------------------------------------------------------------
// Include source density?
bool
pylith::materials::Thermoporoelasticity::useSourceDensity(void) const {
    return _useSourceDensity;
} // useSourceDensity


// ---------------------------------------------------------------------------------------------------------------------
// Include heat source?
void
pylith::materials::Thermoporoelasticity::useHeatSource(const bool value) {
    PYLITH_COMPONENT_DEBUG("useHeatSource(value="<<value<<")");

    _useHeatSource = value;
} // useHeatSource


// ---------------------------------------------------------------------------------------------------------------------
// Include heat source?
bool
pylith::materials::Thermoporoelasticity::useHeatSource(void) const {
    return _useHeatSource;
} // useHeatSource


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain?
void
pylith::materials::Thermoporoelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState(value="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain?
bool
pylith::materials::Thermoporoelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use state variables to update auxiliary fields?
void
pylith::materials::Thermoporoelasticity::useStateVars(const bool value) {
    PYLITH_COMPONENT_DEBUG("useStateVars(value="<<value<<")");

    _useStateVars = value;
} // useStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Use state variables to update auxiliary fields?
bool
pylith::materials::Thermoporoelasticity::useStateVars(void) const {
    return _useStateVars;
} // useStateVars


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Thermoporoelasticity::setBulkRheology(pylith::materials::RheologyThermoporoelasticity* const rheology) {
    PYLITH_COMPONENT_DEBUG("setBulkRheology(rheology="<<rheology<<")");

    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyThermoporoelasticity*
pylith::materials::Thermoporoelasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Thermoporoelasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains required fields
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' subfield in solution for thermoporoelasticity.");
    } // if
    if (!solution.hasSubfield("pressure")) {
        throw std::runtime_error("Cannot find 'pressure' subfield in solution for thermoporoelasticity.");
    } // if
    if (!solution.hasSubfield("temperature")) {
        throw std::runtime_error("Cannot find 'temperature' subfield in solution for thermoporoelasticity.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Thermoporoelasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());
    integrator->createLabelDS(solution, solution.getMesh().getDimension());

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::Thermoporoelasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                              const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    assert(_rheology);

    pylith::materials::AuxiliaryFactoryThermoporoelasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();
    assert(auxiliaryFactory);

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("Thermoporoelasticity auxiliary field");

    assert(_scales);
    auxiliaryFactory->initialize(auxiliaryField, *_scales, domainMesh.getDimension());

    // Add base poroelastic subfields
    auxiliaryFactory->addSolidDensity();
    auxiliaryFactory->addFluidDensity();
    auxiliaryFactory->addFluidViscosity();
    auxiliaryFactory->addPorosity();

    // Add optional fields
    if (_useBodyForce) {
        auxiliaryFactory->addBodyForce();
    } // if
    if (_gravityField) {
        auxiliaryFactory->addGravityField(_gravityField);
    } // if
    if (_useSourceDensity) {
        auxiliaryFactory->addSourceDensity();
    } // if
    if (_useHeatSource) {
        auxiliaryFactory->addHeatSource();
    } // if

    // Add rheology-specific fields
    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::Thermoporoelasticity::createDerivedField(const pylith::topology::Field& solution,
                                                            const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 0) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->setLabel("Thermoporoelasticity derived field");

    assert(_scales);
    _derivedFactory->initialize(derivedField, *_scales, domainMesh.getDimension());
    _derivedFactory->addSubfields();

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    derivedField->allocate();
    derivedField->createOutputVector();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ---------------------------------------------------------------------------------------------------------------------
// Get default PETSc solver options.
pylith::utils::PetscOptions*
pylith::materials::Thermoporoelasticity::getSolverDefaults(const bool isParallel,
                                                           const bool hasFault) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSolverDefaults(isParallel="<<isParallel<<", hasFault="<<hasFault<<")");

    // Use field split for thermoporoelasticity (similar to poroelasticity)
    pylith::utils::PetscOptions* options = new pylith::utils::PetscOptions();

    // TODO: Add appropriate default solver options for thermoporoelasticity
    // This would typically involve a field split approach for the coupled system

    PYLITH_METHOD_RETURN(options);
} // getSolverDefaults


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Thermoporoelasticity::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Thermoporoelasticity::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::Thermoporoelasticity::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::materials::Thermoporoelasticity::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                             const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    const bool hasGravityField = _gravityField != NULL;
    const bool hasTraceStrain = solution.hasSubfield("trace_strain");

    std::vector<ResidualKernels> kernels;

    if (!_useStateVars) {
        // Displacement equation: f1u = stress
        PetscPointFn* f0u = NULL;
        PetscPointFn* f1u = _rheology->getKernelf1u_implicit(coordsys);

        // Pressure equation: f0p = fluid content time derivative, f1p = Darcy flux
        PetscPointFn* f0p = _rheology->getKernelf0p_implicit(coordsys, _useSourceDensity);
        PetscPointFn* f1p = _rheology->getKernelf1p_implicit(coordsys, hasGravityField);

        // Trace strain equation (if present)
        PetscPointFn* f0e = hasTraceStrain ? ThermoporoelasticityKernels::f0e : NULL;
        PetscPointFn* f1e = NULL;

        // Temperature equation: f0T = heat capacity term, f1T = heat flux
        PetscPointFn* f0T = _rheology->getKernelf0T_implicit(coordsys, _useHeatSource);
        PetscPointFn* f1T = _rheology->getKernelf1T_implicit(coordsys);

        if (hasTraceStrain) {
            kernels.resize(4);
            kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
            kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
            kernels[2] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, f0e, f1e);
            kernels[3] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
        } else {
            kernels.resize(3);
            kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
            kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
            kernels[2] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
        }
    } else {
        // State variable formulation
        // Displacement equation
        PetscPointFn* f0u = NULL;
        PetscPointFn* f1u = _rheology->getKernelf1u_implicit(coordsys);

        // Pressure equation
        PetscPointFn* f0p = _rheology->getKernelf0p_implicit(coordsys, _useSourceDensity);
        PetscPointFn* f1p = _rheology->getKernelf1p_implicit(coordsys, hasGravityField);

        // Trace strain equation
        PetscPointFn* f0e = hasTraceStrain ? ThermoporoelasticityKernels::f0e : NULL;
        PetscPointFn* f1e = NULL;

        // Temperature equation
        PetscPointFn* f0T = _rheology->getKernelf0T_implicit(coordsys, _useHeatSource);
        PetscPointFn* f1T = _rheology->getKernelf1T_implicit(coordsys);

        // Velocity equation: f0_v = ∂u/∂t - v = 0
        PetscPointFn* f0v = solution.hasSubfield("velocity") ? ThermoporoelasticityKernels::f0v_implicit : NULL;
        PetscPointFn* f1v = NULL;

        // Pressure_dot equation: f0_pdot = ∂p/∂t - p_dot = 0
        PetscPointFn* f0pdot = solution.hasSubfield("pressure_t") ? ThermoporoelasticityKernels::f0pdot : NULL;
        PetscPointFn* f1pdot = NULL;

        // Trace_strain_dot equation: f0_edot = ∂ε_v/∂t - ε_v_dot = 0
        PetscPointFn* f0edot = solution.hasSubfield("trace_strain_t") ? ThermoporoelasticityKernels::f0edot : NULL;
        PetscPointFn* f1edot = NULL;

        // Temperature_dot equation: f0_Tdot = ∂T/∂t - T_dot = 0
        PetscPointFn* f0Tdot = solution.hasSubfield("temperature_t") ? ThermoporoelasticityKernels::f0Tdot : NULL;
        PetscPointFn* f1Tdot = NULL;

        // Count number of kernels needed
        size_t numKernels = 3; // displacement, pressure, temperature
        if (hasTraceStrain) { numKernels++; }
        if (solution.hasSubfield("velocity")) { numKernels++; }
        if (solution.hasSubfield("pressure_t")) { numKernels++; }
        if (solution.hasSubfield("trace_strain_t")) { numKernels++; }
        if (solution.hasSubfield("temperature_t")) { numKernels++; }

        kernels.resize(numKernels);
        size_t idx = 0;
        kernels[idx++] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        kernels[idx++] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
        if (hasTraceStrain) {
            kernels[idx++] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, f0e, f1e);
        }
        kernels[idx++] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
        if (solution.hasSubfield("velocity")) {
            kernels[idx++] = ResidualKernels("velocity", pylith::feassemble::Integrator::LHS, f0v, f1v);
        }
        if (solution.hasSubfield("pressure_t")) {
            kernels[idx++] = ResidualKernels("pressure_t", pylith::feassemble::Integrator::LHS, f0pdot, f1pdot);
        }
        if (solution.hasSubfield("trace_strain_t")) {
            kernels[idx++] = ResidualKernels("trace_strain_t", pylith::feassemble::Integrator::LHS, f0edot, f1edot);
        }
        if (solution.hasSubfield("temperature_t")) {
            kernels[idx++] = ResidualKernels("temperature_t", pylith::feassemble::Integrator::LHS, f0Tdot, f1Tdot);
        }
    } // if/else _useStateVars

    // Add any MMS body force kernels
    kernels.insert(kernels.end(), _mmsBodyForceKernels.begin(), _mmsBodyForceKernels.end());

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::Thermoporoelasticity::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                             const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    const bool hasTraceStrain = solution.hasSubfield("trace_strain");
    const pylith::feassemble::Integrator::EquationPart equationPart = pylith::feassemble::Integrator::LHS;

    integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

    std::vector<JacobianKernels> kernels;

    if (!_useStateVars) {
        // Get all kernel functions from rheology
        PetscPointJacFn* Jf3uu = _rheology->getKernelJf3uu(coordsys);
        PetscPointJacFn* Jf2up = _rheology->getKernelJf2up(coordsys);
        PetscPointJacFn* Jf2uT = _rheology->getKernelJf2uT(coordsys);
        PetscPointJacFn* Jf0pp = _rheology->getKernelJf0pp(coordsys);
        PetscPointJacFn* Jf3pp = _rheology->getKernelJf3pp(coordsys);
        PetscPointJacFn* Jf0pe = hasTraceStrain ? _rheology->getKernelJf0pe(coordsys) : NULL;
        PetscPointJacFn* Jf0pT = _rheology->getKernelJf0pT(coordsys);
        PetscPointJacFn* Jf0TT = _rheology->getKernelJf0TT(coordsys);
        PetscPointJacFn* Jf3TT = _rheology->getKernelJf3TT(coordsys);

        if (hasTraceStrain) {
            // With trace_strain: 9 kernels
            kernels.resize(9);
            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, NULL, NULL, NULL, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure", equationPart, NULL, NULL, Jf2up, NULL);
            kernels[2] = JacobianKernels("displacement", "temperature", equationPart, NULL, NULL, Jf2uT, NULL);
            kernels[3] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, NULL, NULL, Jf3pp);
            kernels[4] = JacobianKernels("pressure", "trace_strain", equationPart, Jf0pe, NULL, NULL, NULL);
            kernels[5] = JacobianKernels("pressure", "temperature", equationPart, Jf0pT, NULL, NULL, NULL);
            kernels[6] = JacobianKernels("trace_strain", "displacement", equationPart, NULL, ThermoporoelasticityKernels::Jf1eu, NULL, NULL);
            kernels[7] = JacobianKernels("trace_strain", "trace_strain", equationPart, ThermoporoelasticityKernels::Jf0ee, NULL, NULL, NULL);
            kernels[8] = JacobianKernels("temperature", "temperature", equationPart, Jf0TT, NULL, NULL, Jf3TT);
        } else {
            // Without trace_strain: 6 kernels
            kernels.resize(6);
            kernels[0] = JacobianKernels("displacement", "displacement", equationPart, NULL, NULL, NULL, Jf3uu);
            kernels[1] = JacobianKernels("displacement", "pressure", equationPart, NULL, NULL, Jf2up, NULL);
            kernels[2] = JacobianKernels("displacement", "temperature", equationPart, NULL, NULL, Jf2uT, NULL);
            kernels[3] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, NULL, NULL, Jf3pp);
            kernels[4] = JacobianKernels("pressure", "temperature", equationPart, Jf0pT, NULL, NULL, NULL);
            kernels[5] = JacobianKernels("temperature", "temperature", equationPart, Jf0TT, NULL, NULL, Jf3TT);
        }
    } else {
        // State variable formulation
        // Get all kernel functions from rheology
        PetscPointJacFn* Jf3uu = _rheology->getKernelJf3uu(coordsys);
        PetscPointJacFn* Jf2up = _rheology->getKernelJf2up(coordsys);
        PetscPointJacFn* Jf2uT = _rheology->getKernelJf2uT(coordsys);
        PetscPointJacFn* Jf0pp = _rheology->getKernelJf0pp(coordsys);
        PetscPointJacFn* Jf3pp = _rheology->getKernelJf3pp(coordsys);
        PetscPointJacFn* Jf0pe = hasTraceStrain ? _rheology->getKernelJf0pe(coordsys) : NULL;
        PetscPointJacFn* Jf0pT = _rheology->getKernelJf0pT(coordsys);
        PetscPointJacFn* Jf0TT = _rheology->getKernelJf0TT(coordsys);
        PetscPointJacFn* Jf3TT = _rheology->getKernelJf3TT(coordsys);

        // Count number of kernels needed
        size_t numKernels = 6; // base: uu, up, uT, pp, pT, TT
        if (hasTraceStrain) { numKernels += 3; } // pe, eu, ee
        if (solution.hasSubfield("velocity")) { numKernels += 2; } // vu, vv
        if (solution.hasSubfield("pressure_t")) { numKernels += 2; } // pdotp, pdotpdot
        if (solution.hasSubfield("trace_strain_t")) { numKernels += 2; } // edote, edotedot
        if (solution.hasSubfield("temperature_t")) { numKernels += 2; } // TdotT, TdotTdot

        kernels.resize(numKernels);
        size_t idx = 0;

        // Base kernels
        kernels[idx++] = JacobianKernels("displacement", "displacement", equationPart, NULL, NULL, NULL, Jf3uu);
        kernels[idx++] = JacobianKernels("displacement", "pressure", equationPart, NULL, NULL, Jf2up, NULL);
        kernels[idx++] = JacobianKernels("displacement", "temperature", equationPart, NULL, NULL, Jf2uT, NULL);
        kernels[idx++] = JacobianKernels("pressure", "pressure", equationPart, Jf0pp, NULL, NULL, Jf3pp);

        if (hasTraceStrain) {
            kernels[idx++] = JacobianKernels("pressure", "trace_strain", equationPart, Jf0pe, NULL, NULL, NULL);
        }

        kernels[idx++] = JacobianKernels("pressure", "temperature", equationPart, Jf0pT, NULL, NULL, NULL);

        if (hasTraceStrain) {
            kernels[idx++] = JacobianKernels("trace_strain", "displacement", equationPart, NULL, ThermoporoelasticityKernels::Jf1eu, NULL, NULL);
            kernels[idx++] = JacobianKernels("trace_strain", "trace_strain", equationPart, ThermoporoelasticityKernels::Jf0ee, NULL, NULL, NULL);
        }

        kernels[idx++] = JacobianKernels("temperature", "temperature", equationPart, Jf0TT, NULL, NULL, Jf3TT);

        // State variable kernels
        if (solution.hasSubfield("velocity")) {
            kernels[idx++] = JacobianKernels("velocity", "displacement", equationPart, ThermoporoelasticityKernels::Jf0vu, NULL, NULL, NULL);
            kernels[idx++] = JacobianKernels("velocity", "velocity", equationPart, ThermoporoelasticityKernels::Jf0vv, NULL, NULL, NULL);
        }
        if (solution.hasSubfield("pressure_t")) {
            kernels[idx++] = JacobianKernels("pressure_t", "pressure", equationPart, ThermoporoelasticityKernels::Jf0pdotp, NULL, NULL, NULL);
            kernels[idx++] = JacobianKernels("pressure_t", "pressure_t", equationPart, ThermoporoelasticityKernels::Jf0pdotpdot, NULL, NULL, NULL);
        }
        if (solution.hasSubfield("trace_strain_t")) {
            kernels[idx++] = JacobianKernels("trace_strain_t", "trace_strain", equationPart, ThermoporoelasticityKernels::Jf0edote, NULL, NULL, NULL);
            kernels[idx++] = JacobianKernels("trace_strain_t", "trace_strain_t", equationPart, ThermoporoelasticityKernels::Jf0edotedot, NULL, NULL, NULL);
        }
        if (solution.hasSubfield("temperature_t")) {
            kernels[idx++] = JacobianKernels("temperature_t", "temperature", equationPart, ThermoporoelasticityKernels::Jf0TdotT, NULL, NULL, NULL);
            kernels[idx++] = JacobianKernels("temperature_t", "temperature_t", equationPart, ThermoporoelasticityKernels::Jf0TdotTdot, NULL, NULL, NULL);
        }
    } // if/else _useStateVars

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Thermoporoelasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                                 const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 0) {
        PYLITH_METHOD_END;
    } // if

    // Set kernels for derived fields.
    // DerivedFactoryPoroelasticity supports: cauchy_stress, cauchy_strain, bulk_density, water_content
    // We must provide kernels for ALL derived subfields since the kernelsArray is indexed by subfield index.
    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;
    PetscPointFn* bulkDensityKernel = pylith::fekernels::Poroelasticity::bulkDensity_asScalar;

    std::vector<ProjectKernels> kernels(4);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelCauchyStressVector(coordsys));
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);
    kernels[2] = ProjectKernels("bulk_density", bulkDensityKernel);
    kernels[3] = ProjectKernels("water_content", _rheology->getKernelFluidContent(coordsys));

    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
