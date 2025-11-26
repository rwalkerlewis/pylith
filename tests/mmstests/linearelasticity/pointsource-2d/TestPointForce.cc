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

#include "TestPointForce.hh" // Implementation of test class

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/MeshIOPetsc.hh" // USES MeshIOPetsc
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <cassert> // USES assert()


// ------------------------------------------------------------------------------------------------
void
pylith::TestPointForce::setUp(void) {
    MMSTest::setUp();
    _data = NULL;
} // setUp


// ------------------------------------------------------------------------------------------------
void
pylith::TestPointForce::tearDown(void) {
    delete _data;_data = NULL;

    MMSTest::tearDown();
} // tearDown


// ------------------------------------------------------------------------------------------------
void
pylith::TestPointForce::setData(TestPointForce_Data* data) {
    assert(data);
    _data = data;
} // setData


// ------------------------------------------------------------------------------------------------
void
pylith::TestPointForce::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    pythia::journal::debug_t debug("TestPointForce");
    if (debug.state()) {
        debug << pythia::journal::at(__HERE__)
              << "Setting up MMS test for point force "
              << _data->journalName << pythia::journal::endl;
    } // if

    // Set exact solution functions
    assert(_data->exactSolnFns);
    assert(_data->exactSolnDotFns);
    _setExactSolution(_data->exactSolnFns, _data->exactSolnDotFns);

    // Setup problem
    assert(!_problem);
    _problem = new pylith::problems::TimeDependent();assert(_problem);
    _problem->setIdentifier("timedependent");

    assert(_data->numSolnSubfields > 0);
    assert(_data->solnDiscretizations);
    pylith::topology::Field::Discretization* solnDiscretizations =
        const_cast<pylith::topology::Field::Discretization*>(_data->solnDiscretizations);
    _problem->setSolutionDiscretizations(solnDiscretizations, _data->numSolnSubfields);

    _problem->setScales(_data->scales);
    _problem->setFormulation(_data->formulation);
    _isJacobianLinear = _data->isJacobianLinear;

    // Setup mesh
    pylith::topology::Mesh* mesh = NULL;
    if (_data->useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.filename(_data->meshFilename);
        mesh = new pylith::topology::Mesh();
        iohandler.read(mesh);
    } else {
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.filename(_data->meshFilename);
        mesh = new pylith::topology::Mesh();
        iohandler.read(mesh);
    }
    assert(mesh);

    mesh->setCoordSys(&_data->cs);
    pylith::topology::MeshOps::nondimensionalize(mesh, _data->scales);
    _mesh.reset(mesh);

    // Setup material
    _data->material.setFormulation(_data->formulation);
    _data->material.setAuxiliaryFieldDB(&_data->auxDB);
    for (int i = 0; i < _data->numAuxSubfields; ++i) {
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[i];
        _data->material.setAuxiliarySubfieldDiscretization(_data->auxSubfields[i],
                                                           discretization.basisOrder, discretization.quadOrder,
                                                           discretization.dimension,
                                                           discretization.cellBasis, discretization.feSpace,
                                                           discretization.isBasisContinuous);
    } // for
    _data->material.setIdentifier("elasticity");
    _data->material.setName("material-id=24");
    _data->material.setLabelValue(24);
    _data->material.setBulkRheology(&_data->rheology);

    // Setup point source
    _data->pointSource.setFormulation(_data->formulation);
    _data->pointSource.setIdentifier("pointsource");
    _data->pointSource.setName("point-source");
    _data->pointSource.setLabelValue(24); // Same label as material

    std::vector<pylith::materials::Material*> materials(1);
    materials[0] = &_data->material;
    _problem->setMaterials(materials.data(), 1);

    // Add point source to problem (would need integration in Problem class)
    // For now, we test the source in isolation

    // Setup boundary conditions
    for (size_t i = 0; i < _data->bcs.size(); ++i) {
        _data->bcs[i]->setFormulation(_data->formulation);
    }
    _problem->setBoundaryConditions(_data->bcs.data(), _data->bcs.size());

    _problem->preinitialize(*mesh);
    _problem->verifyConfiguration();
    _problem->initialize();

    _t = _data->t;
    _dt = _data->dt;
    _tolerance = _data->tolerance;

    PYLITH_METHOD_END;
} // _initialize


// End of file
