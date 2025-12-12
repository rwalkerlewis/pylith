// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestHeat;
        class TestHeat_Data;
    } // materials
} // pylith

class pylith::materials::TestHeat : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestHeat);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testUseHeatSource);
    CPPUNIT_TEST(testBulkRheology);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test useHeatSource().
    void testUseHeatSource(void);

    /// Test setBulkRheology() and getBulkRheology().
    void testBulkRheology(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestHeat

// End of file
