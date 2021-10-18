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

/**
 * @file tests/libtests/materials/TestAuxiliaryFactoryPoroelasticity.hh
 *
 * @brief C++ TestAuxiliaryFactoryPoroelasticity object.
 *
 * C++ unit testing for AuxiliaryFactoryPoroelasticity.
 */

#if !defined(pylith_materials_testauxiliaryfactoryporoelasticity_hh)
#define pylith_materials_testauxiliaryfactoryporoelasticity_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // HOLDSA AuxiliaryFactoryPoroelasticity
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

#include <map> // USES std::map

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryPoroelasticity;
        class TestAuxiliaryFactoryPoroelasticity_Data;
    } // materials
} // pylith

class pylith::materials::TestAuxiliaryFactoryPoroelasticity : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestAuxiliaryFactoryPoroelasticity);

    CPPUNIT_TEST(testAdd);
    CPPUNIT_TEST(testSetValuesFromDB);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test adding density, body force, and gravity subfields.
    void testAdd(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, auxiliary field, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    AuxiliaryFactoryPoroelasticity* _factory; ///< Test subject.
    TestAuxiliaryFactoryPoroelasticity_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for test subject.

    size_t _auxDim; ///< Topological dimension of auxiliary field.
  
}; // class TestAuxiliaryFactoryPoroelasticity

// =====================================================================================================================
class pylith::materials::TestAuxiliaryFactoryPoroelasticity_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAuxiliaryFactoryPoroelasticity_Data(void);

    /// Destructor
    ~TestAuxiliaryFactoryPoroelasticity_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    size_t dimension; ///< Spatial dimension.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* auxiliaryDB; ///< Spatial database with values for solution.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field spatial database.

}; // class TestAuxiliaryFactoryPoroelasticity_Data

#endif // pylith_materials_testauxiliaryfactoryporoelasticity_hh

// End of file
