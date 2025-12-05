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

#include "TestThermoelasticity.hh" // Implementation of class methods

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/MeshIOPetsc.hh" // USES MeshIOPetsc
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::zeroInitialGuess()
#include "pylith/problems/Physics.hh" // USES Physics

#include "petscts.h" // USES PetscTS

#include "catch2/catch_test_macros.hpp"

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::TestThermoelasticity::TestThermoelasticity(TestThermoelasticity_Data* data) :
    _data(data) {
    _initialize();
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::TestThermoelasticity::~TestThermoelasticity(void) {
    delete _data;_data = NULL;
} // destructor


// ------------------------------------------------------------------------------------------------
// Initialize objects for test.
void
pylith::TestThermoelasticity::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    // Scales for nondimensionalization
    _data->scales.setLengthScale(1.0e+3);      // 1 km
    _data->scales.setTimeScale(1.0e+3);        // 1000 s
    _data->scales.setPressureScale(1.0e+9);    // 1 GPa
    _data->scales.setTemperatureScale(1.0e+3); // 1000 K
    _data->scales.computeDensityScale();

    // Create mesh
    pylith::topology::Mesh* mesh = new pylith::topology::Mesh;
    if (_data->useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.read(mesh);
    } else {
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.setFilename(_data->meshFilename);
        iohandler.setOptions(_data->meshOptions);
        iohandler.read(mesh);
    } // if/else
    mesh->setCoordSys(&_data->cs);
    pylith::topology::MeshOps::nondimensionalize(mesh, _data->scales);
    _mesh = mesh;

    // Setup material
    _data->material.setBulkRheology(&_data->rheology);

    // Set up problem
    pylith::problems::TimeDependent* problem = new pylith::problems::TimeDependent();
    problem->setNormalizer(_data->scales);
    problem->setFormulation(_data->formulation);
    problem->setStartTime(0.0);
    problem->setEndTime(2.0);
    problem->setInitialTimeStep(_data->dt);
    problem->setMaxTimeSteps(1);
    problem->setSolverType(pylith::problems::TimeDependent::NONLINEAR);
    problem->setGravityField(NULL);
    _problem = problem;

    assert(_data->numSolnSubfields > 0);
    assert(_data->solnDiscretizations);
    _problem->defaults().set(pylith::problems::SolutionFactory::displacement(),
                             _data->solnDiscretizations[0].basisOrder,
                             _data->solnDiscretizations[0].quadOrder, true);
    _problem->defaults().set(pylith::problems::SolutionFactory::temperature(),
                             _data->solnDiscretizations[1].basisOrder,
                             _data->solnDiscretizations[1].quadOrder, true);

    // Boundary conditions
    for (size_t i = 0; i < _data->bcs.size(); ++i) {
        _problem->setBoundaryCondition(_data->bcs[i]->getName(), _data->bcs[i]);
    } // for

    // Create solution field (displacement + temperature)
    static const char* subfieldNames[2] = {"displacement", "temperature"};
    pylith::topology::Field* solution = pylith::testing::MMSTest::createSolutionField(*mesh, subfieldNames, 2, _data->solnDiscretizations);
    assert(solution);

    // Set auxiliary data
    _data->material.setLabelName("material-id");
    _data->material.setLabelValue(1);

    assert(_data->numAuxSubfields > 0);
    assert(_data->auxSubfields);
    assert(_data->auxDiscretizations);
    _data->auxDB.addValue("density", _data->auxDB.queryFn("density"), _data->auxDB.units("density"));
    _data->auxDB.addValue("specific_heat", _data->auxDB.queryFn("specific_heat"), _data->auxDB.units("specific_heat"));
    _data->auxDB.addValue("thermal_conductivity", _data->auxDB.queryFn("thermal_conductivity"), _data->auxDB.units("thermal_conductivity"));
    _data->auxDB.addValue("reference_temperature", _data->auxDB.queryFn("reference_temperature"), _data->auxDB.units("reference_temperature"));
    _data->auxDB.addValue("thermal_expansion_coefficient", _data->auxDB.queryFn("thermal_expansion_coefficient"), _data->auxDB.units("thermal_expansion_coefficient"));
    _data->auxDB.addValue("vs", _data->auxDB.queryFn("vs"), _data->auxDB.units("vs"));
    _data->auxDB.addValue("vp", _data->auxDB.queryFn("vp"), _data->auxDB.units("vp"));

    _data->auxDB.setCoordSys(_data->cs);
    _data->material.setAuxiliaryFieldDB(&_data->auxDB);

    for (size_t i = 0; i < _data->numAuxSubfields; ++i) {
        _data->material.setAuxiliarySubfieldDiscretization(_data->auxSubfields[i],
                                                           _data->auxDiscretizations[i].basisOrder,
                                                           _data->auxDiscretizations[i].quadOrder,
                                                           _data->cs.getSpaceDim(),
                                                           _data->auxDiscretizations[i].isBasisContinuous,
                                                           _data->auxDiscretizations[i].feSpace);
    } // for

    _problem->setMaterial("material", &_data->material);

    _problem->preinitialize(*mesh);
    _problem->verifyConfiguration();
    _problem->initialize();
    _solution = solution;

    // Set exact solution
    _setExactSolution();

    // Set test parameters
    _isJacobianLinear = _data->isJacobianLinear;
    _allowZeroResidual = _data->allowZeroResidual;
    _jacobianConvergenceRate = _data->jacobianConvergenceRate;
    _tolerance = _data->tolerance;

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
// Set exact solution and time derivative of solution in domain.
void
pylith::TestThermoelasticity::_setExactSolution(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    const PylithReal t = _data->t;
    PetscErrorCode err;

    // Set exact solution in domain
    err = DMSetAuxiliaryVec(_solution->getDM(), NULL, 0, 0, _exactSolutionVec);REQUIRE(!err);
    err = PetscObjectCompose((PetscObject)_solution->getDM(), "A", (PetscObject)_exactSolutionVec);REQUIRE(!err);

    // Set exact solution for displacement field
    PetscVec exactDisp = NULL;
    err = VecDuplicate(_solution->getLocalVector(), &exactDisp);REQUIRE(!err);
    err = PetscObjectSetName((PetscObject)exactDisp, "exact_displacement");REQUIRE(!err);

    if (_data->exactSolnFns) {
        err = DMProjectFunction(_solution->getDM(), t, _data->exactSolnFns, NULL, INSERT_VALUES, exactDisp);REQUIRE(!err);
    } // if
    err = VecCopy(exactDisp, _exactSolutionVec);REQUIRE(!err);
    err = VecDestroy(&exactDisp);REQUIRE(!err);

    // Set exact solution time derivative
    if (_data->exactSolnDotFns) {
        PetscVec exactDot = NULL;
        err = VecDuplicate(_solution->getLocalVector(), &exactDot);REQUIRE(!err);
        err = PetscObjectSetName((PetscObject)exactDot, "exact_solution_dot");REQUIRE(!err);
        err = DMProjectFunction(_solution->getDM(), t, _data->exactSolnDotFns, NULL, INSERT_VALUES, exactDot);REQUIRE(!err);
        err = VecCopy(exactDot, _exactSolutionDotVec);REQUIRE(!err);
        err = VecDestroy(&exactDot);REQUIRE(!err);
    } // if

    PYLITH_METHOD_END;
} // _setExactSolution


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::TestThermoelasticity_Data::TestThermoelasticity_Data(void) :
    journalName("TestThermoelasticity"),
    spaceDim(2),
    meshFilename(NULL),
    meshOptions(""),
    boundaryLabel("boundary"),
    useAsciiMesh(true),
    jacobianConvergenceRate(1.0),
    tolerance(1.0e-9),
    isJacobianLinear(true),
    allowZeroResidual(true),
    t(0.0),
    dt(0.05),
    formulation(pylith::problems::Physics::QUASISTATIC),
    numSolnSubfields(0),
    solnDiscretizations(NULL),
    exactSolnFns(NULL),
    exactSolnDotFns(NULL),
    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::TestThermoelasticity_Data::~TestThermoelasticity_Data(void) {}


// End of file
