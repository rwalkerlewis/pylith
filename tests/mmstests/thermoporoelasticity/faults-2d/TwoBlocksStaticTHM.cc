// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file tests/mmstests/thermoporoelasticity/faults-2d/TwoBlocksStaticTHM.cc
 *
 * Square domain of sides 12.0 km with one through-going fault running
 * through the domain at y=+2km. The two opposing sides each
 * move as rigid blocks with 1.5 m of right-lateral slip.
 * Pressure and temperature are uniform at reference values (no thermal or pore pressure strain).
 */

#include <portinfo>

#include "TwoBlocksStaticTHM.hh" // Implementation of cases

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin
#include "pylith/faults/KinSrcStep.hh" // USES KinSrcStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Thermoporoelasticity.hh" // USES Thermoporoelasticity
#include "pylith/materials/IsotropicLinearThermoporoelasticity.hh" // USES IsotropicLinearThermoporoelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Mesh.hh" // USES pylith::topology::Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "pylith/scales/Scales.hh" // USES Scales
#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

namespace pylith {
    class _TwoBlocksStaticTHM;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_TwoBlocksStaticTHM {
    static pylith::scales::Scales scales;
    static const double AMPLITUDE; // nondimensional
    static const double X_FAULT; // dimensional
    static const double REFERENCE_TEMP; // dimensional
    static const double REFERENCE_PRES; // dimensional (zero - no pore pressure effect)

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
        return -AMPLITUDE * scales.getLengthScale();
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

    // Pressure (constant at zero - no pore pressure effect)
    static double fluid_pressure(const double x,
                                 const double y) {
        return 0.0;
    } // fluid_pressure

    // Trace strain (zero for rigid block motion)
    static double trace_strain(const double x,
                               const double y) {
        return 0.0;
    } // trace_strain

    // Temperature (constant at reference temperature)
    static double temperature(const double x,
                              const double y) {
        const PylithReal temperatureScale = scales.getTemperatureScale();
        return REFERENCE_TEMP / temperatureScale;
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

        return PETSC_SUCCESS;
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

        return PETSC_SUCCESS;
    } // solnkernel_trace_strain

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
    TestFaultKinThermoporoelasticity_Data* createData(void) {
        TestFaultKinThermoporoelasticity_Data* data = new TestFaultKinThermoporoelasticity_Data();assert(data);

        data->journalName = "TwoBlocksStaticTHM";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.

        scales = data->scales;

        // solnDiscretizations set in derived class.

        data->matNumAuxSubfields = 14;
        static const char* _matAuxSubfields[14] = {
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
            "specific_heat"
        };
        data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[14] = {
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
        data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        data->matAuxDB.addValue("solid_density", solid_density, density_units());
        data->matAuxDB.addValue("fluid_density", fluid_density, density_units());
        data->matAuxDB.addValue("fluid_viscosity", fluid_viscosity, viscosity_units());
        data->matAuxDB.addValue("porosity", porosity, porosity_units());
        data->matAuxDB.addValue("shear_modulus", shear_modulus, modulus_units());
        data->matAuxDB.addValue("drained_bulk_modulus", drained_bulk_modulus, modulus_units());
        data->matAuxDB.addValue("biot_coefficient", biot_coefficient, biot_coefficient_units());
        data->matAuxDB.addValue("fluid_bulk_modulus", fluid_bulk_modulus, modulus_units());
        data->matAuxDB.addValue("isotropic_permeability", isotropic_permeability, permeability_units());
        data->matAuxDB.addValue("reference_temperature", reference_temperature, temperature_units());
        data->matAuxDB.addValue("thermal_expansion_coefficient", thermal_expansion_coefficient, thermal_expansion_units());
        data->matAuxDB.addValue("fluid_thermal_expansion", fluid_thermal_expansion, thermal_expansion_units());
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
            pylith::materials::Thermoporoelasticity* material = new pylith::materials::Thermoporoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoporoelasticity");
            material->setName("material-id=10");
            material->setLabelValue(10);
            material->setBulkRheology(&data->rheology);
            data->materials[0] = material;
        } // xneg
        { // mid
            pylith::materials::Thermoporoelasticity* material = new pylith::materials::Thermoporoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoporoelasticity");
            material->setName("material-id=20");
            material->setLabelValue(20);
            material->setBulkRheology(&data->rheology);
            data->materials[1] = material;
        } // mid
        { // xpos
            pylith::materials::Thermoporoelasticity* material = new pylith::materials::Thermoporoelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setIdentifier("thermoporoelasticity");
            material->setName("material-id=15");
            material->setLabelValue(15);
            material->setBulkRheology(&data->rheology);
            data->materials[2] = material;
        } // xpos

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt constrainedScalar[1] = {0};
        static const PylithInt numConstrainedDOF = 2;
        static const PylithInt numConstrainedScalar = 1;
        data->bcs.resize(8);
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
        { // boundary_xpos (pressure)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[2] = bc;
        } // boundary_xpos
        { // boundary_xneg (pressure)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[3] = bc;
        } // boundary_xneg
        { // boundary_xpos (temperature)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[4] = bc;
        } // boundary_xpos
        { // boundary_xneg (temperature)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("temperature");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_temperature);
            data->bcs[5] = bc;
        } // boundary_xneg
        { // boundary_yneg (pressure - no flow)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[6] = bc;
        } // boundary_yneg
        { // boundary_ypos (pressure - no flow)
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedScalar, numConstrainedScalar);
            bc->setUserFn(solnkernel_fluid_pressure);
            data->bcs[7] = bc;
        } // boundary_ypos

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

        data->numSolnSubfieldsDomain = 4;
        data->numSolnSubfieldsFault = 1;
        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[5] = {
            solnkernel_disp,
            solnkernel_fluid_pressure,
            solnkernel_trace_strain,
            solnkernel_temperature,
            solnkernel_lagrangemultiplier,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // _TwoBlocksStaticTHM
pylith::scales::Scales pylith::_TwoBlocksStaticTHM::scales;
const double pylith::_TwoBlocksStaticTHM::AMPLITUDE = 3.0;
const double pylith::_TwoBlocksStaticTHM::X_FAULT = +2.0e+3;
const double pylith::_TwoBlocksStaticTHM::REFERENCE_TEMP = 300.0;
const double pylith::_TwoBlocksStaticTHM::REFERENCE_PRES = 0.0;

// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoporoelasticity_Data*
pylith::TwoBlocksStaticTHM::TriP2P1P1P1(void) {
    TestFaultKinThermoporoelasticity_Data* data = pylith::_TwoBlocksStaticTHM::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    data->allowZeroResidual = true;

    assert(4 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
        pylith::topology::Field::Discretization(1, 2), // trace_strain
        pylith::topology::Field::Discretization(1, 2), // temperature
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[14] = {
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
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    return data;
} // TriP2P1P1P1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoporoelasticity_Data*
pylith::TwoBlocksStaticTHM::TriP3P2P2P2(void) {
    TestFaultKinThermoporoelasticity_Data* data = pylith::_TwoBlocksStaticTHM::createData();assert(data);

    data->meshFilename = "data/tri.mesh";
    data->allowZeroResidual = true;

    assert(4 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
        pylith::topology::Field::Discretization(2, 3), // trace_strain
        pylith::topology::Field::Discretization(2, 3), // temperature
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[14] = {
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
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    return data;
} // TriP3P2P2P2


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoporoelasticity_Data*
pylith::TwoBlocksStaticTHM::QuadQ2Q1Q1Q1(void) {
    TestFaultKinThermoporoelasticity_Data* data = pylith::_TwoBlocksStaticTHM::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    data->allowZeroResidual = true;

    assert(4 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // pressure
        pylith::topology::Field::Discretization(1, 2), // trace_strain
        pylith::topology::Field::Discretization(1, 2), // temperature
        pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[14] = {
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
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 2), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    return data;
} // QuadQ2Q1Q1Q1


// ------------------------------------------------------------------------------------------------
pylith::TestFaultKinThermoporoelasticity_Data*
pylith::TwoBlocksStaticTHM::QuadQ3Q2Q2Q2(void) {
    TestFaultKinThermoporoelasticity_Data* data = pylith::_TwoBlocksStaticTHM::createData();assert(data);

    data->meshFilename = "data/quad.mesh";
    data->allowZeroResidual = true;

    assert(4 == data->numSolnSubfieldsDomain);
    assert(1 == data->numSolnSubfieldsFault);
    static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // pressure
        pylith::topology::Field::Discretization(2, 3), // trace_strain
        pylith::topology::Field::Discretization(2, 3), // temperature
        pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _matAuxDiscretizations[14] = {
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
    data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

    static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        pylith::topology::Field::Discretization(0, 3), // slip
    };
    data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

    return data;
} // QuadQ3Q2Q2Q2


// End of file
