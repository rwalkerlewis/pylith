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

#include "UniformHeat2D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

// ------------------------------------------------------------------------------------------------
namespace pylith {
    class _UniformHeat2D;
}
class pylith::_UniformHeat2D {
private:

    // Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Specific heat capacity
    static double specific_heat(const double x,
                                const double y) {
        return 900.0;
    } // specific_heat

    static const char* specific_heat_units(void) {
        return "J/(kg*K)";
    } // specific_heat_units

    // Thermal conductivity
    static double thermal_conductivity(const double x,
                                       const double y) {
        return 3.0;
    } // thermal_conductivity

    static const char* thermal_conductivity_units(void) {
        return "W/(m*K)";
    } // thermal_conductivity_units

    // Solution subfields.

    // Temperature gradient components
    static double grad_T_x(void) {
        return 10.0; // K/m
    } // grad_T_x

    static double grad_T_y(void) {
        return 20.0; // K/m
    } // grad_T_y

    // Temperature: T = T0 + a*x + b*y (linear temperature field)
    // For this MMS test, we use T = 300 + 10*x + 20*y
    // This gives a constant gradient, so the Laplacian is zero
    static double temperature(const double x,
                              const double y) {
        return 300.0 + grad_T_x()*x + grad_T_y()*y;
    } // temperature

    static const char* temperature_units(void) {
        return "K";
    } // temperature_units

    static
    PetscErrorCode solnkernel_temperature(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(1 == numComponents);
        assert(s);

        s[0] = temperature(x[0], x[1]);

        return PETSC_SUCCESS;
    } // solnkernel_temperature

public:

    static
    TestHeat_Data* createData(void) {
        TestHeat_Data* data = new TestHeat_Data();assert(data);

        data->journalName = "UniformHeat2D";

        data->isJacobianLinear = true;
        data->allowZeroResidual = true; // Linear temperature field -> Laplacian = 0

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "specific_heat", "thermal_conductivity"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // specific_heat
            pylith::topology::Field::Discretization(0, 1), // thermal_conductivity
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("specific_heat", specific_heat, specific_heat_units());
        data->auxDB.addValue("thermal_conductivity", thermal_conductivity, thermal_conductivity_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->material.useHeatSource(false);

        data->material.setIdentifier("heat");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[1] = {0};
        static const PylithInt numConstrained = 1;
        data->bcs.resize(1);
        pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("temperature");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_temperature);
        data->bcs[0] = bc;

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[1] = {
            solnkernel_temperature,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // _UniformHeat2D

// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::TriP1(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::TriP2(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // specific_heat
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::TriP3(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // specific_heat
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::QuadQ1(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::QuadQ2(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // specific_heat
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestHeat_Data*
pylith::UniformHeat2D::QuadQ3(void) {
    TestHeat_Data* data = pylith::_UniformHeat2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // specific_heat
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// End of file
