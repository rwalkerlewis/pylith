// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file tests/mmstests/thermoelasticity/faults-2d/TwoBlocksStaticThermal.cc
 *
 * Square domain of sides 12.0 km with one through-going fault running
 * through the domain at y=+2km. The two opposing sides each
 * move as rigid blocks with 1.5 m of right-lateral slip.
 * Temperature is uniform at reference temperature (no thermal strain).
 */

#include <portinfo>

#include "TwoBlocksStaticThermal.hh" // Implementation of cases

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcStep.hh" // USES KinSrcStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Thermoelasticity.hh" // USES Thermoelasticity
#include "pylith/materials/IsotropicLinearThermoelasticity.hh" // USES IsotropicLinearThermoelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Mesh.hh" // USES pylith::topology::Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales

namespace pylith {
    class _TwoBlocksStaticThermal;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_TwoBlocksStaticThermal {
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double RIGIDITY_SCALE;
    static const double TEMPERATURE_SCALE;
    static const double AMPLITUDE; // nondimensional
    static const double X_FAULT; // nondimensional
    static const double REFERENCE_TEMP; // dimensional

    // Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y) {
        return 3000.0;
    } // vs

    static const char* vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y) {
        return sqrt(3.0)*vs(x,y);
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

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

    static double thermal_conductivity(const double x,
                                       const double y) {
        return 3.0;
    } // thermal_conductivity

    static const char* thermal_conductivity_units(void) {
        return "W/(m*K)";
    } // thermal_conductivity_units

    static double specific_heat(const double x,
                                const double y) {
        return 800.0;
    } // specific_heat

    static const char* specific_heat_units(void) {
        return "J/(kg*K)";
    } // specific_heat_units

    // Kinematic rupture auxiliary components.

    // Initiation time
    static double initiation_time(const double x,
                                  const double y) {
        return 0.0;
    } // initiation_time

    static const char* time_units(void) {
        return "s";
    } // time_units

    // Slip
    static double finalslip_opening(const double x,
                                    const double y) {
        return 0.0;
    } // slip_opening

    static double finalslip_leftlateral(const double x,
                                        const double y) {
        return -AMPLITUDE * LENGTH_SCALE;
    } // slip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         PetscInt flag) {
        const double disp = 0.5*AMPLITUDE;
        if (!flag) {
            return x < X_FAULT ? +disp : -disp;
        } else {
            return flag < 0 ? +disp : -disp;
        } // if/else
    } // disp_y

    // Temperature (constant at reference temperature)
    static double temperature(const double x,
                              const double y) {
        return REFERENCE_TEMP / TEMPERATURE_SCALE;
    } // temperature

    static double faulttraction_x(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_x

    static double faulttraction_y(const double x,
                                  const double y) {
        return 0.0;
    } // faulttraction_y

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1]);
        PetscInt flag = 0;
        if (context) {
            PetscErrorCode err = PETSC_SUCCESS;
            PetscDM dmMesh = PetscDM(context);
            PetscInt cell = 0;
            err = DMPlexGetActivePoint(dmMesh, &cell);PYLITH_CHECK_ERROR(err);

            double centroid[3] = {0.0, 0.0, 0.0};
            err = DMPlexComputeCellGeometryFVM(dmMesh, cell, NULL, centroid, NULL);PYLITH_CHECK_ERROR(err);
            flag = centroid[0] < X_FAULT ? -1 : +1;
        } // if
        s[1] = disp_y(x[0], x[1], flag);

        return PETSC_SUCCESS;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_temperature(PetscInt spaceDim,
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

    static PetscErrorCode solnkernel_lagrangemultiplier(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = faulttraction_x(x[0], x[1]);
        s[1] = faulttraction_y(x[0], x[1]);

        return PETSC_SUCCESS;
    } // solnkernel_lagrangemultiplier

public:

    static
    TestFaultKinThermoelasticity_Data* createData(void) {
        TestFaultKinThermoelasticity_Data* data = new TestFaultKinThermoelasticity_Data();assert(data);

        data->journalName = "TwoBlocksStaticThermal";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.

        data->scales.setLengthScale(LENGTH_SCALE);
        data->scales.setTimeScale(TIME_SCALE);
        data->scales.setRigidityScale(RIGIDITY_SCALE);
        data->scales.setTemperatureScale(TEMPERATURE_SCALE);
        data->scales.computeDensityScale();

        // solnDiscretizations set in derived class.

        data->matNumAuxSubfields = 7;
        static const char* _matAuxSubfields[7] = {
            "density",
            "shear_modulus",
            "bulk_modulus",
            "reference_temperature",
            "thermal_expansion_coefficient",
            "thermal_conductivity",
            "specific_heat"
        };
        data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // reference_temperature
            pylith::topology::Field::Discretization(0, 1), // thermal_expansion_coefficient
            pylith::topology::Field::Discretization(0, 1), // thermal_conductivity
            pylith::topology::Field::Discretization(0, 1), // specific_heat
        };
        data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        data->matAuxDB.addValue("density", density, density_units());
        data->matAuxDB.addValue("vp", vp, vp_units());
        data->matAuxDB.addValue("vs", vs, vs_units());
        data->matAuxDB.addValue("reference_temperature", reference_temperature, temperature_units());
        data->matAuxDB.addValue("thermal_expansion_coefficient", thermal_expansion_coefficient, thermal_expansion_units());
        data->matAuxDB.addValue("thermal_conductivity", thermal_conductivity, thermal_conductivity_units());
        data->matAuxDB.addValue("specific_heat", specific_heat, specific_heat_units());
        data->matAuxDB.setCoordSys(data->cs);

        assert(!data->kinSrc);
        data->kinSrc = new pylith::faults::KinSrcStep();assert(data->kinSrc);
        data->kinSrc->setOriginTime(0.0);
        data->faultAuxDB.addValue("initiation_time", initiation_time, time_units());
        data->faultAuxDB.addValue("final_slip_opening", finalslip_opening, slip_units());
        data->faultAuxDB.addValue("final_slip_left_lateral", finalslip_leftlateral, slip_units());
        data->faultAuxDB.setCoordSys(data->cs);

        data->faultNumAuxSubfields = 1;
        static const char* _faultAuxSubfields[1] = { "slip" };
        data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 1), // slip
        };
        data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Materials
        data->materials.resize(3);
        { // xneg
            pylith::materials::Thermoelasticity* material = new pylith::materials::Thermoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoelasticity");
            material->setName("material-id=10");
            material->setLabelValue(10);
            material->setBulkRheology(&data->rheology);
            data->materials[0] = material;
        } // xneg
        { // mid
            pylith::materials::Thermoelasticity* material = new pylith::materials::Thermoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoelasticity");
            material->setName("material-id=20");
            material->setLabelValue(20);
            material->setBulkRheology(&data->rheology);
            data->materials[1] = material;
        } // mid
        { // xpos
            pylith::materials::Thermoelasticity* material = new pylith::materials::Thermoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoelasticity");
            material->setName("material-id=15");
            material->setLabelValue(15);
            material->setBulkRheology(&data->rheology);
            data->materials[2] = material;
        } // xpos

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt constrainedScalar[1] = {0};
        static const PylithInt numConstrainedDOF = 2;
        static const PylithInt numConstrainedScalar = 1;
        data->bcs.resize(4);
        { // boundary_xpos (displacement)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrainedDOF);
            bc->setUserFn(solnkernel_disp);
            data->bcs[0] = bc;
        } // boundary_xpos
        { // boundary_xneg (displacement)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrainedDOF);
            bc->setUserFn(solnkernel_disp);
            data->bcs[1] = bc;
        } // boundary_xneg
        { // boundary_xpos (temperature)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[2] = bc;
        } // boundary_xpos
        { // boundary_xneg (temperature)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[3] = bc;
        } // boundary_xneg

        // Faults
        data->faults.resize(1);
        { // xpos
            pylith::faults::FaultCohesiveKin* fault = new pylith::faults::FaultCohesiveKin();
            fault->setCohesiveLabelValue(100);
            fault->setSurfaceLabelName("fault_xpos_faces");

            const int numRuptures = 1;
            const char* ruptureNames[1] = { "rupture" };
            pylith::faults::KinSrc* ruptures[1] = { data->kinSrc };
            fault->setEqRuptures(ruptureNames, numRuptures, ruptures, numRuptures);
            data->faults[0] = fault;
        } // xpos

        pylith::utils::PetscOptions options;
        options.add("-fieldsplit_displacement_pc_type", "lu");
        options.override ();

        data->numSolnSubfieldsDomain = 2;
        data->numSolnSubfieldsFault = 1;
        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[3] = {
            solnkernel_disp,
            solnkernel_temperature,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // _TwoBlocksStaticThermal
const double pylith::_TwoBlocksStaticThermal::LENGTH_SCALE = 1.0;
const double pylith::_TwoBlocksStaticThermal::TIME_SCALE = 2.0;
const double pylith::_TwoBlocksStaticThermal::RIGIDITY_SCALE = 2.0e+6;
const double pylith::_TwoBlocksStaticThermal::TEMPERATURE_SCALE = 1.0e+3;
const double pylith::_TwoBlocksStaticThermal::AMPLITUDE = 3.0;
const double pylith::_TwoBlocksStaticThermal::X_FAULT = +2.0e+3;
const double pylith::_TwoBlocksStaticThermal::REFERENCE_TEMP = 300.0;

// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::TriP1(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    data->allowZeroResidual = true;

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // temperature
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::TriP2(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // temperature
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::TriP3(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // temperature
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::QuadQ1(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // temperature
        pylith::topology::Field::Discretization(1, 1, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::QuadQ2(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // specific_heat
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // temperature
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoelasticity_Data*
pylith::TwoBlocksStaticThermal::QuadQ3(void) {
    TestFaultKinThermoelasticity_Data* data = pylith::_TwoBlocksStaticThermal::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        pylith::topology::Field::Discretization(0, 3), // reference_temperature
        pylith::topology::Field::Discretization(0, 3), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 3), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 3), // specific_heat
    };
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    assert(2 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(3, 3), // temperature
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ3


// End of file
