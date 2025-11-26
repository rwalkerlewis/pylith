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

#include "pylith/testing/MMSTest.hh" // ISA MMSTest

#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity

#include "pylith/sources/PointForce.hh" // USES PointForce

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/problems/Physics.hh" // USES FormulationEnum
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    class TestPointForce;
    class TestPointForce_Data;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::TestPointForce_Data {
public:

    const char* journalName; ///< Name for journal.
    bool isJacobianLinear; ///< True if Jacobian is linear.

    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* boundaryLabel; ///< Label for boundary.
    bool useAsciiMesh; ///< True if using ASCII mesh; false for Gmsh.

    pylith::scales::Scales scales; ///< Scales for nondimensionalization.

    pylith::problems::Physics::FormulationEnum formulation; ///< Formulation

    PylithReal t; ///< Time at which to check solution.
    PylithReal dt; ///< Time step.
    PylithReal tolerance; ///< Tolerance for test.

    int numSolnSubfields; ///< Number of solution subfields.
    pylith::topology::Field::Discretization const* solnDiscretizations; ///< Discretization for solution subfields.

    int numAuxSubfields; ///< Number of material auxiliary subfields.
    const char* const* auxSubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization const* auxDiscretizations; ///< Discretization for auxiliary subfields.

    spatialdata::spatialdb::UserFunctionDB auxDB; ///< Spatial database with auxiliary field.
    spatialdata::geocoords::CSCart cs; ///< Coordinate system.

    pylith::materials::Elasticity material; ///< Elasticity material.
    pylith::materials::IsotropicLinearElasticity rheology; ///< Linear elastic rheology.

    pylith::sources::PointForce pointSource; ///< Point force source.

    std::vector<pylith::bc::BoundaryCondition*> bcs; ///< Boundary conditions.

    pylith::testing::MMSTest::solution_fn* exactSolnFns; ///< Exact solution functions.
    pylith::testing::MMSTest::solution_fn* exactSolnDotFns; ///< Time derivatives of exact solution functions.

}; // TestPointForce_Data


// ------------------------------------------------------------------------------------------------
class pylith::TestPointForce : public pylith::testing::MMSTest {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Set data.
    void setData(TestPointForce_Data* data);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestPointForce_Data* _data; ///< Data for testing.

}; // TestPointForce

// End of file
