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

#include "TestHeat.hh" // Implementation of class methods

#include "pylith/materials/Heat.hh" // USES Heat
#include "pylith/materials/IsotropicHeat.hh" // USES IsotropicHeat

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestHeat);

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestHeat::setUp(void) {
    _mesh = new pylith::topology::Mesh();
    CPPUNIT_ASSERT(_mesh);
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestHeat::tearDown(void) {
    delete _mesh;
    _mesh = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::materials::TestHeat::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    Heat heat;

    CPPUNIT_ASSERT_EQUAL(false, heat.useHeatSource());
    CPPUNIT_ASSERT(!heat.getBulkRheology());

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test useHeatSource().
void
pylith::materials::TestHeat::testUseHeatSource(void) {
    PYLITH_METHOD_BEGIN;

    Heat heat;

    CPPUNIT_ASSERT_EQUAL(false, heat.useHeatSource());

    heat.useHeatSource(true);
    CPPUNIT_ASSERT_EQUAL(true, heat.useHeatSource());

    heat.useHeatSource(false);
    CPPUNIT_ASSERT_EQUAL(false, heat.useHeatSource());

    PYLITH_METHOD_END;
} // testUseHeatSource


// ----------------------------------------------------------------------
// Test setBulkRheology() and getBulkRheology().
void
pylith::materials::TestHeat::testBulkRheology(void) {
    PYLITH_METHOD_BEGIN;

    Heat heat;
    IsotropicHeat rheology;

    CPPUNIT_ASSERT(!heat.getBulkRheology());

    heat.setBulkRheology(&rheology);
    CPPUNIT_ASSERT_EQUAL(&rheology, heat.getBulkRheology());

    PYLITH_METHOD_END;
} // testBulkRheology


// End of file
