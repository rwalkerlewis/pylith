// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/source/WellboreSource.hh" // implementation of object methods

#include "pylith/source/AuxiliaryFactoryWellboreSource.hh" // USES AuxiliaryFactoryWellboreSource
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/WellboreSource.hh" // USES WellboreSource kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::source::WellboreSource::WellboreSource(void) :
    _useInertia(false) {
    pylith::utils::PyreComponent::setName("wellboresource");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::source::WellboreSource::~WellboreSource(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::source::WellboreSource::deallocate(void) {
    Material::deallocate();
} // deallocate



// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::source::WellboreSource::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getLabel()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("pressure")) {
        throw std::runtime_error("Cannot find 'pressure' field in solution; required for 'WellboreSource'.");
    } // if
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_END;
} // verifyConfiguration


// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::source::WellboreSource::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getLabel()<<")");

    PetscDM dm = solution.dmMesh();assert(dm);
    // transform points of source to mesh coordinates in python
    // DM from solution
    PetscSF sfPoints;
    Vec vecPoints;
    err = VecCreateMPIWithArray(PetscObjectComm((PetscObject) dm), ,_numPoints, _numPoints, _pointCoords, vecPoints);PYLITH_CHECK_ERROR(err);
    err = DMLocatePoints(dm, vecPoints, DM_POINTLOCATION_NONE, sfPoints);PYLITH_CHECK_ERROR(err);


    pylith::feassemble:IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setLabelName(pylith::topology::Mesh::getCellsLabelName());
    integrator->setLabelValue(getMaterialId());

    switch (_formulation) {
    case QUASISTATIC:
        _useInertia = false;
        break;
    case DYNAMIC:
        _useInertia = true;
        break;
    case DYNAMIC_IMEX:
        _useInertia = true;
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    _setKernelsLHSResidual(integrator, solution);
    _setKernelsLHSJacobian(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::source::WellboreSource::createAuxiliaryField(const pylith::topology::Field& solution,
                                                    const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getLabel()<<", domainMesh="<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->setLabel("WellboreSource auxiliary field");

    assert(_rheology);
    pylith::source::AuxiliaryFactoryWellboreSource* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    // add in aux specific to peaceman
    auxiliaryFactory->addFluidDensity(); // 0
    auxiliaryFactory->addFluidViscosity(); // 1
    auxiliaryFactory->addIsotropicPermeability(); // 2
    auxiliaryFactory->addWellboreRadius(); // 3
    auxiliaryFactory->addWellboreLength(); // 4
    auxiliaryFactory->addWellborePressure(); // 5
    auxiliaryFactory->addWellboreCharacter(); // 6
    auxiliaryFactory->addElementDimensions(); // 7

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
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::source::WellboreSource::_getAuxiliaryFactory(void) {
    return _auxiliaryFactory();
} // _getAuxiliaryFactory

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::source::WellboreSource::_setKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                      const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSResidual(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();

    std::vector<ResidualKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {

        // Pressure
        const PetscPointFunc f0p = pylith::fekernels::WellboreSource::f0p;
        const PetscPointFunc f1p = NULL;        

        kernels.resize(1);
        kernels[0] = ResidualKernels("pressure", f0p, f1p);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX: {
        break;
    } // DYNAMIC
    case DYNAMIC: {
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setKernelsLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::source::WellboreSource::_setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                      const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsLHSJacobian(integrator="<<integrator<<", solution="<<solution.getLabel()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().getCoordSys();

    std::vector<JacobianKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointJac Jf0pp = pylith::fekernels::WellboreSource::Jf0pp;
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        kernels.resize(1);
        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        break;
    } // QUASISTATIC
    case DYNAMIC:
    case DYNAMIC_IMEX: {
        break;
    } // DYNAMIC_IMEX
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsLHSJacobian(kernels);

    PYLITH_METHOD_END;
} // _setKernelsLHSJacobian


// End of file
