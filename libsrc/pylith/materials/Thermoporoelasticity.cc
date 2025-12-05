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
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::fekernels::Thermoporoelasticity ThermoporoelasticityKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Thermoporoelasticity::Thermoporoelasticity(void) :
    _useBodyForce(false),
    _useSourceDensity(false),
    _useHeatSource(false),
    _useReferenceState(false),
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

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

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

    assert(_normalizer);
    _derivedFactory->initialize(derivedField, *_normalizer, domainMesh.getDimension());
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

    std::vector<ResidualKernels> kernels;

    // Displacement equation
    // f1u: stress
    const PetscPointFn f1u = _rheology->getKernelf1u_implicit(coordsys);
    kernels.resize(1);
    kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, NULL, f1u);
    integrator->setKernelsResidual(kernels);

    // Pressure equation
    // f0p: fluid content time derivative
    // f1p: Darcy flux
    const PetscPointFn f0p = _rheology->getKernelf0p_implicit(coordsys, _useSourceDensity);
    const PetscPointFn f1p = _rheology->getKernelf1p_implicit(coordsys, hasGravityField);
    kernels.resize(1);
    kernels[0] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, f1p);
    integrator->setKernelsResidual(kernels);

    // Trace strain equation (if present)
    if (solution.hasSubfield("trace_strain")) {
        kernels.resize(1);
        kernels[0] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, 
                                     ThermoporoelasticityKernels::f0e, NULL);
        integrator->setKernelsResidual(kernels);
    } // if

    // Temperature equation
    // f0T: heat capacity term
    // f1T: heat flux
    const PetscPointFn f0T = _rheology->getKernelf0T_implicit(coordsys, _useHeatSource);
    const PetscPointFn f1T = _rheology->getKernelf1T_implicit(coordsys);
    kernels.resize(1);
    kernels[0] = ResidualKernels("temperature", pylith::feassemble::Integrator::LHS, f0T, f1T);
    integrator->setKernelsResidual(kernels);

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

    std::vector<JacobianKernels> kernels;

    // Displacement-displacement (elastic stiffness)
    const PetscPointJacFn Jf3uu = _rheology->getKernelJf3uu(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("displacement", "displacement", pylith::feassemble::Integrator::LHS,
                                 NULL, NULL, NULL, Jf3uu);
    integrator->setKernelsJacobian(kernels);

    // Displacement-pressure (Biot coupling)
    const PetscPointJacFn Jf2up = _rheology->getKernelJf2up(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("displacement", "pressure", pylith::feassemble::Integrator::LHS,
                                 NULL, NULL, Jf2up, NULL);
    integrator->setKernelsJacobian(kernels);

    // Displacement-temperature (thermal coupling)
    const PetscPointJacFn Jf2uT = _rheology->getKernelJf2uT(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("displacement", "temperature", pylith::feassemble::Integrator::LHS,
                                 NULL, NULL, Jf2uT, NULL);
    integrator->setKernelsJacobian(kernels);

    // Pressure-pressure (storage + Darcy)
    const PetscPointJacFn Jf0pp = _rheology->getKernelJf0pp(coordsys);
    const PetscPointJacFn Jf3pp = _rheology->getKernelJf3pp(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("pressure", "pressure", pylith::feassemble::Integrator::LHS,
                                 Jf0pp, NULL, NULL, Jf3pp);
    integrator->setKernelsJacobian(kernels);

    // Pressure-trace_strain (Biot coupling)
    if (solution.hasSubfield("trace_strain")) {
        const PetscPointJacFn Jf0pe = _rheology->getKernelJf0pe(coordsys);
        kernels.resize(1);
        kernels[0] = JacobianKernels("pressure", "trace_strain", pylith::feassemble::Integrator::LHS,
                                     Jf0pe, NULL, NULL, NULL);
        integrator->setKernelsJacobian(kernels);
    } // if

    // Pressure-temperature (thermal coupling in fluid content)
    const PetscPointJacFn Jf0pT = _rheology->getKernelJf0pT(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("pressure", "temperature", pylith::feassemble::Integrator::LHS,
                                 Jf0pT, NULL, NULL, NULL);
    integrator->setKernelsJacobian(kernels);

    // Trace_strain equation Jacobians (if present)
    if (solution.hasSubfield("trace_strain")) {
        // trace_strain - displacement
        kernels.resize(1);
        kernels[0] = JacobianKernels("trace_strain", "displacement", pylith::feassemble::Integrator::LHS,
                                     NULL, ThermoporoelasticityKernels::Jf1eu, NULL, NULL);
        integrator->setKernelsJacobian(kernels);

        // trace_strain - trace_strain
        kernels.resize(1);
        kernels[0] = JacobianKernels("trace_strain", "trace_strain", pylith::feassemble::Integrator::LHS,
                                     ThermoporoelasticityKernels::Jf0ee, NULL, NULL, NULL);
        integrator->setKernelsJacobian(kernels);
    } // if

    // Temperature-temperature (heat capacity + conductivity)
    const PetscPointJacFn Jf0TT = _rheology->getKernelJf0TT(coordsys);
    const PetscPointJacFn Jf3TT = _rheology->getKernelJf3TT(coordsys);
    kernels.resize(1);
    kernels[0] = JacobianKernels("temperature", "temperature", pylith::feassemble::Integrator::LHS,
                                 Jf0TT, NULL, NULL, Jf3TT);
    integrator->setKernelsJacobian(kernels);

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

    std::vector<ProjectKernels> kernels(3);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelCauchyStressVector(coordsys));
    kernels[1] = ProjectKernels("fluid_content", _rheology->getKernelFluidContent(coordsys));
    kernels[2] = ProjectKernels("heat_flux", _rheology->getKernelHeatFluxVector(coordsys));

    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
