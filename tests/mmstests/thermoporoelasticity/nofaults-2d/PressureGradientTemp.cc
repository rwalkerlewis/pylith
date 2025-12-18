// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file tests/mmstests/thermoporoelasticity/nofaults-2d/PressureGradientTemp.cc
 *
 * MMS test case for thermoporoelasticity with a pressure gradient and
 * uniform temperature at reference temperature.
 *
 * The domain is 8km x 8km centered at the origin.
 * We apply a pressure gradient in the x-direction with uniform temperature
 * at the reference temperature (so thermal strain contribution is zero).
 */

#include <portinfo>

#include "PressureGradientTemp.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

namespace pylith {
    class _PressureGradientTemp;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_PressureGradientTemp {
    static pylith::scales::Scales scales;
    static const double PRESSURE; // dimensional
    static const double X_MAX; // dimensional
    static const double REFERENCE_TEMP; // dimensional

    // Density
    static double solid_density(const double x,
                                const double y) {
        return 2500.0;
    } // solid_density

    static double fluid_density(const double x,
                                const double y) {
        return 1000.0;
    } // fluid_density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Fluid viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 1.0e-3;
    } // fluid_viscosity

    static const char* viscosity_units(void) {
        return "Pa*s";
    } // viscosity_units

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.02;
    } // porosity

    static const char* porosity_units(void) {
        return "none";
    } // porosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 3.0e+10;
    } // shear_modulus

    static const char* modulus_units(void) {
        return "Pa";
    } // modulus_units

    // Drained bulk modulus
    static double drained_bulk_modulus(const double x,
                                       const double y) {
        return 8.0e+10;
    } // drained_bulk_modulus

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 0.8;
    } // biot_coefficient

    static const char* biot_coefficient_units(void) {
        return "none";
    } // biot_coefficient_units

    // Fluid modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y) {
        return 2.0e+9;
    } // fluid_bulk_modulus

    // Permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 1.0e-14;
    } // isotropic_permeability

    static const char* permeability_units(void) {
        return "m*m";
    } // permeability_units

    // Thermal properties
    static double reference_temperature(const double x,
                                        const double y) {
        return REFERENCE_TEMP;
    } // reference_temperature

    static const char* temperature_units(void) {
        return "K";
    } // temperature_units

    static double thermal_expansion_coefficient(const double x,
                                                const double y) {
        return 1.0e-5;
    } // thermal_expansion_coefficient

    static const char* thermal_expansion_units(void) {
        return "1/K";
    } // thermal_expansion_units

    static double fluid_thermal_expansion(const double x,
                                          const double y) {
        return 2.1e-4;
    } // fluid_thermal_expansion

    static double thermal_conductivity(const double x,
                                       const double y) {
        return 3.0;
    } // thermal_conductivity

    static const char* thermal_conductivity_units(void) {
        return "watt/m/K";
    } // thermal_conductivity_units

    static double specific_heat(const double x,
                                const double y) {
        return 800.0;
    } // specific_heat

    static const char* specific_heat_units(void) {
        return "J/(kg*K)";
    } // specific_heat_units

    // Solution subfields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        const PylithReal lengthScale = scales.getLengthScale();
        const PylithReal rigidityScale = scales.getRigidityScale();
        const PylithReal fluidPressureScale = pylith::scales::ElasticityScales::getFluidPressureScale(scales);

        const double muN = shear_modulus(x, y) / rigidityScale;
        const double lambdaN = drained_bulk_modulus(x, y) / rigidityScale - 2.0/3.0 * muN;
        const double alpha = biot_coefficient(x, y);
        return -0.5 * alpha  * (PRESSURE / fluidPressureScale) / (lambdaN + 2.0*muN) * (x*x / (X_MAX / lengthScale));
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y

    // Pressure
    static double fluid_pressure(const double x,
                                 const double y) {
        const PylithReal lengthScale = scales.getLengthScale();
        const PylithReal fluidPressureScale = pylith::scales::ElasticityScales::getFluidPressureScale(scales);

        return (PRESSURE / fluidPressureScale) * (1.0 - x / (X_MAX / lengthScale));
    } // fluid_pressure

    // Trace strain
    static double trace_strain(const double x,
                               const double y) {
        const PylithReal lengthScale = scales.getLengthScale();
        const PylithReal rigidityScale = scales.getRigidityScale();
        const PylithReal fluidPressureScale = pylith::scales::ElasticityScales::getFluidPressureScale(scales);

        const double muN = shear_modulus(x, y) / rigidityScale;
        const double lambdaN = drained_bulk_modulus(x, y) / rigidityScale - 2.0/3.0 * muN;
        const double alpha = biot_coefficient(x, y);
        return -alpha  * (PRESSURE / fluidPressureScale)  / (lambdaN + 2.0*muN) * (x / (X_MAX / lengthScale));
    } // trace_strain

    // Temperature (constant at reference temperature)
    static double temperature(const double x,
                              const double y) {
        const PylithReal temperatureScale = scales.getTemperatureScale();
        return REFERENCE_TEMP / temperatureScale;
    } // temperature

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_fluid_pressure(PetscInt spaceDim,
                                                    PetscReal t,
                                                    const PetscReal x[],
                                                    PetscInt numComponents,
                                                    PetscScalar* s,
                                                    void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = fluid_pressure(x[0], x[1]);

        return 0;
    } // solnkernel_fluid_pressure

    static PetscErrorCode solnkernel_trace_strain(PetscInt spaceDim,
                                                  PetscReal t,
                                                  const PetscReal x[],
                                                  PetscInt numComponents,
                                                  PetscScalar* s,
                                                  void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = trace_strain(x[0], x[1]);

        return 0;
    } // solnkernel_trace_strain

    static PetscErrorCode solnkernel_temperature(PetscInt spaceDim,
                                                 PetscReal t,
                                                 const PetscReal x[],
                                                 PetscInt numComponents,
                                                 PetscScalar* s,
                                                 void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = temperature(x[0], x[1]);

        return 0;
    } // solnkernel_temperature

    static PetscErrorCode solnkernel_velocity(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* s,
                                              void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);

        s[0] = 0.0;
        s[1] = 0.0;

        return 0;
    } // solnkernel_velocity

    static PetscErrorCode solnkernel_fluid_pressure_dot(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_fluid_pressure_dot

    static PetscErrorCode solnkernel_trace_strain_dot(PetscInt spaceDim,
                                                      PetscReal t,
                                                      const PetscReal x[],
                                                      PetscInt numComponents,
                                                      PetscScalar* s,
                                                      void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_trace_strain_dot

    static PetscErrorCode solnkernel_temperature_dot(PetscInt spaceDim,
                                                     PetscReal t,
                                                     const PetscReal x[],
                                                     PetscInt numComponents,
                                                     PetscScalar* s,
                                                     void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(s);

        s[0] = 0.0;

        return 0;
    } // solnkernel_temperature_dot

public:

    static
    TestThermoporoelasticity_Data* createData(void) {
        TestThermoporoelasticity_Data* data = new TestThermoporoelasticity_Data();assert(data);

        data->journalName = "PressureGradientTemp";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        scales = data->scales;

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 14;
        static const char* _auxSubfields[14] = { // order must match order of subfields in auxiliary field
            "solid_density",
            "fluid_density",
            "fluid_viscosity",
            "porosity",
            "shear_modulus",
            "drained_bulk_modulus",
            "biot_coefficient",
            "biot_modulus",
            "isotropic_permeability",
            "reference_temperature",
            "thermal_expansion_coefficient",
            "fluid_thermal_expansion",
            "thermal_conductivity",
            "specific_heat",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
            pylith::topology::Field::Discretization(0, 1), // solid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // biot_coefficient
            pylith::topology::Field::Discretization(0, 1), // biot_modulus
            pylith::topology::Field::Discretization(0, 1), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 1), // reference_temperature
            pylith::topology::Field::Discretization(0, 1), // thermal_expansion_coefficient
            pylith::topology::Field::Discretization(0, 1), // fluid_thermal_expansion
            pylith::topology::Field::Discretization(0, 1), // thermal_conductivity
            pylith::topology::Field::Discretization(0, 1), // specific_heat
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("solid_density", solid_density, density_units());
        data->auxDB.addValue("fluid_density", fluid_density, density_units());
        data->auxDB.addValue("fluid_viscosity", fluid_viscosity, viscosity_units());
        data->auxDB.addValue("porosity", porosity, porosity_units());
        data->auxDB.addValue("shear_modulus", shear_modulus, modulus_units());
        data->auxDB.addValue("drained_bulk_modulus", drained_bulk_modulus, modulus_units());
        data->auxDB.addValue("biot_coefficient", biot_coefficient, modulus_units());
        data->auxDB.addValue("fluid_bulk_modulus", fluid_bulk_modulus, modulus_units());
        data->auxDB.addValue("isotropic_permeability", isotropic_permeability, permeability_units());
        data->auxDB.addValue("reference_temperature", reference_temperature, temperature_units());
        data->auxDB.addValue("thermal_expansion_coefficient", thermal_expansion_coefficient, thermal_expansion_units());
        data->auxDB.addValue("fluid_thermal_expansion", fluid_thermal_expansion, thermal_expansion_units());
        data->auxDB.addValue("thermal_conductivity", thermal_conductivity, thermal_conductivity_units());
        data->auxDB.addValue("specific_heat", specific_heat, specific_heat_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->rheology.useReferenceState(false);

        data->material.setIdentifier("thermoporoelasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        static const PylithInt constrainedX[1] = { 0 };
        static const PylithInt constrainedY[1] = { 1 };
        static const PylithInt constrainedScalar[1] = { 0 };
        static const PylithInt numConstrained = 1;
        data->bcs.resize(8);
        { // Displacement -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        }
        { // Displacement +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedX, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        }
        { // Displacement -y
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[2] = bc;
        }
        { // Displacement +y
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedY, numConstrained);
            bc->setUserFn(solnkernel_disp);
            data->bcs[3] = bc;
        }
        { // Pressure -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrained);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[4] = bc;
        }
        { // Pressure +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrained);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[5] = bc;
        }
        { // Temperature -x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrained);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[6] = bc;
        }
        { // Temperature +x
            pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrained);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[7] = bc;
        }

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[4] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
            solnkernel_temperature,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

    static
    TestThermoporoelasticity_Data* createDataStateVars(void) {
        TestThermoporoelasticity_Data* data = createData();

        data->material.useStateVars(true);

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[8] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
            solnkernel_temperature,
            solnkernel_velocity,
            solnkernel_fluid_pressure_dot,
            solnkernel_trace_strain_dot,
            solnkernel_temperature_dot,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);

        return data;
    } // createDataStateVars

}; // PressureGradientTemp
pylith::scales::Scales pylith::_PressureGradientTemp::scales;
const double pylith::_PressureGradientTemp::PRESSURE = 4.0e+6;
const double pylith::_PressureGradientTemp::X_MAX = 8.0e+3;
const double pylith::_PressureGradientTemp::REFERENCE_TEMP = 300.0;

// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::TriP2P1P1P1(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 4;
    static const pylith::topology::Field::Discretization _solnDiscretizations[4] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(1, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2P1P1P1


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::TriP3P2P2P2(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 4;
    static const pylith::topology::Field::Discretization _solnDiscretizations[4] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(2, 3), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3P2P2P2


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::QuadQ2Q1Q1Q1(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 4;
    static const pylith::topology::Field::Discretization _solnDiscretizations[4] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(1, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1Q1


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::QuadQ3Q2Q2Q2(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 4;
    static const pylith::topology::Field::Discretization _solnDiscretizations[4] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(2, 3), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3Q2Q2Q2


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::TriP2P1P1P1_StateVars(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createDataStateVars();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 8;
    static const pylith::topology::Field::Discretization _solnDiscretizations[8] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(1, 2), // temperature
        pylith::topology::Field::Discretization(2, 2), // velocity
        pylith::topology::Field::Discretization(1, 2), // fluid pressure dot
        pylith::topology::Field::Discretization(1, 2), // trace strain dot
        pylith::topology::Field::Discretization(1, 2), // temperature dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2P1P1P1_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::TriP3P2P2P2_StateVars(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createDataStateVars();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 8;
    static const pylith::topology::Field::Discretization _solnDiscretizations[8] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(2, 3), // temperature
        pylith::topology::Field::Discretization(3, 3), // velocity
        pylith::topology::Field::Discretization(2, 3), // fluid pressure dot
        pylith::topology::Field::Discretization(2, 3), // trace strain dot
        pylith::topology::Field::Discretization(2, 3), // temperature dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3P2P2P2_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createDataStateVars();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 8;
    static const pylith::topology::Field::Discretization _solnDiscretizations[8] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(1, 2), // fluid pressure
        pylith::topology::Field::Discretization(1, 2), // trace strain
        pylith::topology::Field::Discretization(1, 2), // temperature
        pylith::topology::Field::Discretization(2, 2), // velocity
        pylith::topology::Field::Discretization(1, 2), // fluid pressure dot
        pylith::topology::Field::Discretization(1, 2), // trace strain dot
        pylith::topology::Field::Discretization(1, 2), // temperature dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 2), // solid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_density
        pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 2), // porosity
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // biot_coefficient
        pylith::topology::Field::Discretization(0, 2), // biot_modulus
        pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2Q1Q1Q1_StateVars


// ------------------------------------------------------------------------------------------------
pylith::TestThermoporoelasticity_Data*
pylith::PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars(void) {
    TestThermoporoelasticity_Data* data = pylith::_PressureGradientTemp::createDataStateVars();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 8;
    static const pylith::topology::Field::Discretization _solnDiscretizations[8] = {
        pylith::topology::Field::Discretization(3, 3), // displacement
        pylith::topology::Field::Discretization(2, 3), // fluid pressure
        pylith::topology::Field::Discretization(2, 3), // trace strain
        pylith::topology::Field::Discretization(2, 3), // temperature
        pylith::topology::Field::Discretization(3, 3), // velocity
        pylith::topology::Field::Discretization(2, 3), // fluid pressure dot
        pylith::topology::Field::Discretization(2, 3), // trace strain dot
        pylith::topology::Field::Discretization(2, 3), // temperature dot
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[14] = {
        pylith::topology::Field::Discretization(0, 3), // solid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_density
        pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        pylith::topology::Field::Discretization(0, 3), // porosity
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // biot_coefficient
        pylith::topology::Field::Discretization(0, 3), // biot_modulus
        pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // fluid_thermal_expansion
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3Q2Q2Q2_StateVars


// End of file
