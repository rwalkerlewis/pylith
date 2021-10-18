// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestAuxiliaryFactoryPoroelasticity.hh" // Implementation of class methods

#include "pylith/materials/AuxiliaryFactoryPoroelasticity.hh" // Test subject

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/testing/FieldTester.hh" // USES FieldTester
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestAuxiliaryFactoryPoroelasticity::setUp(void) {
    PYLITH_METHOD_BEGIN;
    _data = new TestAuxiliaryFactoryPoroelasticity_Data();CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT(_data->normalizer);
    _data->normalizer->setLengthScale(1.0e+03);
    _data->normalizer->setTimeScale(2.0);
    _data->normalizer->setDensityScale(3.0e+3);
    _data->normalizer->setPressureScale(2.25e+10);

    pylith::topology::Field::SubfieldInfo info;
    pylith::string_vector componentNames;

    // body_force
    componentNames.resize(3);
    componentNames[0] = "body_force_x";
    componentNames[1] = "body_force_y";
    componentNames[2] = "body_force_z";
    info.description = pylith::topology::Field::Description(
        "body_force",
        "body_force",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->normalizer->getPressureScale() / _data->normalizer->getLengthScale()
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _auxDim, _auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, false
        );
    info.index = 1;
    _data->subfields["body_force"] = info;

    // gravity_field
    componentNames.resize(3);
    componentNames[0] = "gravitational_acceleration_x";
    componentNames[1] = "gravitational_acceleration_y";
    componentNames[2] = "gravitational_acceleration_z";
    info.description = pylith::topology::Field::Description(
        "gravitational_acceleration",
        "gravitational_acceleration",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::VECTOR,
        _data->normalizer->getLengthScale() / pow(_data->normalizer->getTimeScale(), 2)
        );
    info.fe = pylith::topology::Field::Discretization(
        2, 2, _auxDim, _auxDim, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 2;
    _data->subfields["gravitational_acceleration"] = info;

    // porosity
    componentNames.resize(1);
    componentNames[0] = "porosity";
    info.description = pylith::topology::Field::Description(
        "porosity",
        "porosity",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        1,
        pylith::topology::FieldQuery::validatorNonnegative
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["porosity"] = info;

    // solid_density
    componentNames.resize(1);
    componentNames[0] = "solid_density";
    info.description = pylith::topology::Field::Description(
        "solid_density",
        "solid_density",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getDensityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["solid_density"] = info;

    // fluid_density
    componentNames.resize(1);
    componentNames[0] = "fluid_density";
    info.description = pylith::topology::Field::Description(
        "fluid_density",
        "fluid_density",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getDensityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["fluid_density"] = info;

    // fluid_viscosity
    componentNames.resize(1);
    componentNames[0] = "fluid_viscosity";
    info.description = pylith::topology::Field::Description(
        "fluid_viscosity",
        "fluid_viscosity",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getPressureScale() * _data->normalizer->getTimeScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["fluid_viscosity"] = info;

    // source_density
    componentNames.resize(1);
    componentNames[0] = "source_density";
    info.description = pylith::topology::Field::Description(
        "source_density",
        "source_density",
        componentNames,
        componentNames.size(),
        pylith::topology::Field::SCALAR,
        _data->normalizer->getDensityScale(),
        pylith::topology::FieldQuery::validatorPositive
        );
    info.fe = pylith::topology::Field::Discretization(
        1, 2, _auxDim, 1, false, pylith::topology::Field::DEFAULT_BASIS, pylith::topology::Field::POLYNOMIAL_SPACE, true
        );
    info.index = 0;
    _data->subfields["source_density"] = info;

    PYLITH_METHOD_END;
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Tear down testing data.
void
pylith::materials::TestAuxiliaryFactoryPoroelasticity::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _factory;_factory = NULL;
    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _auxiliaryField;_auxiliaryField = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ---------------------------------------------------------------------------------------------------------------------
// Test adding density, body force, and gravity subfields.
void
pylith::materials::TestAuxiliaryFactoryPoroelasticity::testAdd(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);
    CPPUNIT_ASSERT(_data);

    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("solid_density"));
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("fluid_density"));
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("fluid_viscosity"));
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("porosity"));    
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("body_force"));
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("gravitational_acceleration"));
    CPPUNIT_ASSERT(!_auxiliaryField->hasSubfield("source_density"));    

    _factory->addSolidDensity();
    _factory->addFluidDensity();
    _factory->addFluidViscosity();
    _factory->addPorosity();    
    _factory->addBodyForce();
    CPPUNIT_ASSERT(_data->gravityField);
    _factory->addGravityField(_data->gravityField);
    _factory->addSourceDensity();

    CPPUNIT_ASSERT(_data->normalizer);

    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["solid_density"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["fluid_density"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["fluid_viscosity"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["porosity"]);    
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["body_force"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["gravitational_acceleration"]);
    pylith::testing::FieldTester::checkSubfieldInfo(*_auxiliaryField, _data->subfields["source_density"]);    

    PYLITH_METHOD_END;
} // testAdd


// ---------------------------------------------------------------------------------------------------------------------
// Test setValues().
void
pylith::materials::TestAuxiliaryFactoryPoroelasticity::testSetValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_factory);

    _factory->addSolidDensity();
    _factory->addFluidDensity();
    _factory->addFluidViscosity();
    _factory->addBodyForce();
    CPPUNIT_ASSERT(_data->gravityField);
    _factory->addGravityField(_data->gravityField);
    _factory->addSourceDensity();
    _auxiliaryField->subfieldsSetup();
    _auxiliaryField->createDiscretization();
    _auxiliaryField->allocate();

    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    _factory->setValuesFromDB();
    pylith::testing::FieldTester::checkFieldWithDB(*_auxiliaryField, _data->auxiliaryDB, _data->normalizer->getLengthScale());

    PYLITH_METHOD_END;
} // testSetValues


// ---------------------------------------------------------------------------------------------------------------------
// Initialze mesh, coordinate system, auxiliary field, and factory.
void
pylith::materials::TestAuxiliaryFactoryPoroelasticity::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    iohandler.read(_mesh);

    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any cells.",
                           pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    CPPUNIT_ASSERT_MESSAGE("Test mesh does not contain any vertices.",
                           pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);

    // Setup coordinates.
    _mesh->setCoordSys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _auxiliaryField = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_auxiliaryField);
    _auxiliaryField->setLabel("auxiliary");

    _factory = new AuxiliaryFactoryPoroelasticity();
    CPPUNIT_ASSERT(_data->auxiliaryDB);
    _factory->setQueryDB(_data->auxiliaryDB);
    typedef std::map<std::string, pylith::topology::Field::SubfieldInfo>::const_iterator subfield_iter;
    for (subfield_iter iter = _data->subfields.begin(); iter != _data->subfields.end(); ++iter) {
        const char* subfieldName = iter->first.c_str();
        const pylith::topology::Field::Discretization& fe = iter->second.fe;
        _factory->setSubfieldDiscretization(subfieldName, fe.basisOrder, fe.quadOrder, fe.dimension, fe.isFaultOnly, fe.cellBasis,
                                            fe.feSpace, fe.isBasisContinuous);
    } // for
    CPPUNIT_ASSERT(_data->normalizer);
    _factory->initialize(_auxiliaryField, *_data->normalizer, _data->dimension);

    PYLITH_METHOD_END;
} // _initialize


// ---------------------------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryPoroelasticity_Data::TestAuxiliaryFactoryPoroelasticity_Data(void) :
    meshFilename(NULL),
    cs(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    auxiliaryDB(new spatialdata::spatialdb::UserFunctionDB),
    gravityField(new spatialdata::spatialdb::GravityField)
{}


// ---------------------------------------------------------------------------------------------------------------------
pylith::materials::TestAuxiliaryFactoryPoroelasticity_Data::~TestAuxiliaryFactoryPoroelasticity_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxiliaryDB;auxiliaryDB = NULL;
    delete gravityField;gravityField = NULL;
} // destructor


// End of file
