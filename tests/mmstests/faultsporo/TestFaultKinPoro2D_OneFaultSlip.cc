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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file tests/mmstests/faults/TestFaultKinPoro2D_OneFaultSlip.cc
 *
 * Square domain of sides 12.0 km with a single through-going fault at y=+0 km.
 * Prescribed slip of two meters left lateral.
 */

#include <portinfo>

#include "TestFaultKinPoro.hh" // ISA TestFaultKinPoro

#include "pylith/faults/FaultCohesiveKinPoro.hh" // USES FaultCohesiveKinPoro
#include "pylith/faults/KinSrcPoroStep.hh" // USES KinSrcPoroStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn
#include "pylith/bc/NeumannUserFn.hh" // USES NeumannUserFn

#include "pylith/topology/Mesh.hh" // USES pylith::topology::Mesh::cells_label_name
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestFaultKinPoro2D_OneFaultSlip;

        class TestFaultKinPoro2D_OneFaultSlip_TriP1;
        class TestFaultKinPoro2D_OneFaultSlip_TriP2;
        class TestFaultKinPoro2D_OneFaultSlip_TriP3;
        class TestFaultKinPoro2D_OneFaultSlip_TriP4;

        class TestFaultKinPoro2D_OneFaultSlip_QuadQ1;
        class TestFaultKinPoro2D_OneFaultSlip_QuadQ2;
        class TestFaultKinPoro2D_OneFaultSlip_QuadQ3;
        class TestFaultKinPoro2D_OneFaultSlip_QuadQ4;
    } // tests/mmstests
} // pylith

// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip :
    public pylith::mmstests::TestFaultKinPoro {
    // Spatial database user functions for material auxiiliary subfields (includes derived fields).

    // Fluid density
    static double fluid_density(const double x,
                                const double y) {
        return 2500.0;
    } // fluid_density

    static const char* fluid_density_units(void) {
        return "kg/m**3";
    } // fluid_density_units

    // Solid density
    static double solid_density(const double x,
                                const double y) {
        return 2500.0;
    } // solid_density

    static const char* solid_density_units(void) {
        return "kg/m**3";
    } // solid_density_units

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.1;
    } // porosity

    static const char* porosity_units(void) {
        return "None";
    } // porosity_units

    // Fluid viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 1.0;
    } // fluid_viscosity

    static const char* fluid_viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 0.5;
    } // shear_modulus

    static const char* shear_modulus_units(void) {
        return "Pa";
    } // shear_modulus_units

    // Drained bulk modulus
    static double drained_bulk_modulus(const double x,
                                       const double y) {
        return 4.0/3.0;
    } // drained_bulk_modulus

    static const char* drained_bulk_modulus_units(void) {
        return "Pa";
    } // drained_bulk_modulus_units

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 0.5;
    } // biot_coefficient

    static const char* biot_coefficient_units(void) {
        return "None";
    } // biot_coefficient_units

    // Biot modulus
    static double biot_modulus(const double x,
                               const double y) {
        return 1.0;
    } // biot_coefficient

    static const char* biot_modulus_units(void) {
        return "Pa**-1";
    } // biot_modulus_units

    // Fluid bulk modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y) {
        return 0.5;
    } // fluid_bulk_modulus

    static const char* fluid_bulk_modulus_units(void) {
        return "Pa";
    } // fluid_bulk_modulus_units

    // Solid bulk modulus
    static double solid_bulk_modulus(const double x,
                                     const double y) {
        return 0.5;
    } // solid_bulk_modulus

    static const char* solid_bulk_modulus_units(void) {
        return "Pa";
    } // solid_bulk_modulus_units

    // Isotropic permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 1.0;
    } // isotropic_permeability

    static const char* isotropic_permeability_units(void) {
        return "m**2";
    } // isotropic_permeability_units

    // Spatial database user functions for fault auxiiliary subfields (includes derived fields).

    // Fault Thickness
    static double thickness(const double x,
                            const double y) {
        return 1.0;
    } // thickness

    static const char* thickness_units(void) {
        return "m";
    } // thickness_units

    // Fault Porosity
    static double fault_porosity(const double x,
                                 const double y) {
        return 0.1;
    } // fault_porosity

    static const char* fault_porosity_units(void) {
        return "None";
    } // fault_porosity_units

    // Beta p
    static double beta_p(const double x,
                         const double y) {
        return 1.0;
    } // beta_p

    static const char* beta_p_units(void) {
        return "Pa**-1";
    } // beta_p_units

    // Beta sigma
    static double beta_sigma(const double x,
                             const double y) {
        return 1.0;
    } // beta_sigma

    static const char* beta_sigma_units(void) {
        return "Pa**-1";
    } // beta_sigma_units

    // Fault permeability
    static double fault_permeability_xx(const double x,
                                        const double y) {
        return 1.0;
    } // fault_permeability_xx

    static double fault_permeability_yy(const double x,
                                        const double y) {
        return 1.0;
    } // fault_permeability_yy

    static double fault_permeability_zz(const double x,
                                        const double y) {
        return 1.0;
    } // fault_permeability_zz

    static double fault_permeability_xy(const double x,
                                        const double y) {
        return 0.0;
    } // fault_permeability_xy

    static const char* fault_permeability_units(void) {
        return "m**2";
    } // fault_permeability_units

    // Fault fluid viscosity
    static double fault_fluid_viscosity(const double x,
                                        const double y) {
        return 1.0;
    } // fluid_viscosity

    static const char* fault_fluid_viscosity_units(void) {
        return "Pa*s";
    } // fault_fluid_viscosity_units

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
        return 2.0;
    } // slip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Displacement
    static double disp_x_neg(const double x,
                         const double y,
                         const double t) {
        return x * x - 1;
    } // disp_x_neg
    static double disp_x_pos(const double x,
                         const double y,
                         const double t) {
        return 1 - x * x;
    } // disp_x_pos

    static double disp_y_neg(const double x,
                         const double y,
                         const double t) {
        return x * x;
    } // disp_y_neg
    static double disp_y_pos(const double x,
                         const double y,
                         const double t) {
        return -x * x;
    } // disp_y_pos

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static double pressure(const double x,
                           const double y,
                           const double t) {
        return t * (x * x - 1);
    }

    static const char* pressure_units(void) {
        return "Pa";
    } // pressure_units

    static double trace_strain(const double x,
                               const double y,
                               const double t) {
        return (x > 0) ? -2 * x : 2 * x;

    } // trace_strain

    static const char* trace_strain_units(void) {
        return "None";
    } // trace_strain_units

    static double lagrange_multiplier_x(const double x,
                                        const double y,
                                        const double t) {
        return 0;
    } // lagrange_multiplier_x

    static double lagrange_multiplier_y(const double x,
                                        const double y,
                                        const double t) {
        return 0.5 * t;
    } // lagrange_multiplier_y

    static const char* lagrange_multiplier_units(void) {
        return "None";
    } // lagrange_multiplier_units

    static double fault_pressure(const double x,
                                 const double y,
                                 const double t) {
        return 5.5 * t + 0.25 * x * x;
    } // fault_pressure

    static const char* fault_pressure(void) {
        return "Pa";
    } // fault_pressure

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        // This is a HACK. We should make separate DSes for the two sides of the fault, which take separate exact solutions
        DM dm = (DM) context;
        PetscInt point;
        DMPlexGetActivePoint(dm, &point);
        if (point < 12) {
          s[0] = disp_x_neg(x[0], x[1], t);
          s[1] = disp_y_neg(x[0], x[1], t);
        } else {
          s[0] = disp_x_pos(x[0], x[1], t);
          s[1] = disp_y_pos(x[0], x[1], t);
        }

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_pres(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = pressure(x[0], x[1], t);

        return 0;
    } // solnkernel_pres

    static PetscErrorCode solnkernel_tracestrain(PetscInt spaceDim,
                                                 PetscReal t,
                                                 const PetscReal x[],
                                                 PetscInt numComponents,
                                                 PetscScalar* s,
                                                 void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = trace_strain(x[0], x[1], t);

        return 0;
    } // solnkernel_tracestrain

    static PetscErrorCode solnkernel_lagrangemultiplier(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt numComponents,
                                                        PetscScalar* s,
                                                        void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = lagrange_multiplier_x(x[0], x[1], t);
        s[1] = lagrange_multiplier_y(x[0], x[1], t);

        return 0;
    } // solnkernel_lagrangemultiplier

    static PetscErrorCode solnkernel_faultpressure(PetscInt spaceDim,
                                                   PetscReal t,
                                                   const PetscReal x[],
                                                   PetscInt numComponents,
                                                   PetscScalar* s,
                                                   void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = fault_pressure(x[0], x[1], t);

        return 0;
    } // solnkernel_faultpressure

    // ----------------------------------------------------------------------
    // f0 function for displacement forcing.
    static void f0u(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
                    const PylithScalar s[],
                    const PylithScalar s_t[],
                    const PylithScalar s_x[],
                    const PylithInt aOff[],
                    const PylithInt aOff_x[],
                    const PylithScalar a[],
                    const PylithScalar a_t[],
                    const PylithScalar a_x[],
                    const PylithReal t,
                    const PylithScalar x[],
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        f0[0] = (x[0] > 0) ? 4.0 + x[0] * t : -4.0 + x[0] * t;
        f0[1] = (x[0] > 0) ? 1.0 : -1.0;

    } // f0u

    // ----------------------------------------------------------------------
    // f0 function for pressure forcing.
    static void f0p(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
                    const PylithScalar s[],
                    const PylithScalar s_t[],
                    const PylithScalar s_x[],
                    const PylithInt aOff[],
                    const PylithInt aOff_x[],
                    const PylithScalar a[],
                    const PylithScalar a_t[],
                    const PylithScalar a_x[],
                    const PylithReal t,
                    const PylithScalar x[],
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        f0[0] = 1.0 * (x[0] * x[0] - 1) - 1.9 * t;
    } // f0p

    // ----------------------------------------------------------------------
    // f0 function for trace strain forcing.
    static void f0e(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
                    const PylithScalar s[],
                    const PylithScalar s_t[],
                    const PylithScalar s_x[],
                    const PylithInt aOff[],
                    const PylithInt aOff_x[],
                    const PylithScalar a[],
                    const PylithScalar a_t[],
                    const PylithScalar a_x[],
                    const PylithReal t,
                    const PylithScalar x[],
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        f0[0] = 0.0;
    } // f0e

protected:

    void setUp(void) {
        TestFaultKinPoro::setUp();

        // Overwrite component name for control of journals at test level.
        GenericComponent::setName("TestFaultKinPoro2D_OneFaultSlip");

        CPPUNIT_ASSERT(!_data);
        _data = new TestFaultKinPoro_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0);
        _data->normalizer->setTimeScale(1.0);
        _data->normalizer->setPressureScale(1.0);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 5.0;
        _data->timeStep = 1.0;

        // solnDiscretizations set in derived class.

        _data->matNumAuxSubfields = 9;
        static const char *_matAuxSubfields[9] = {"solid_density", "fluid_density", "fluid_viscosity", "porosity", "shear_modulus", "drained_bulk_modulus", "biot_coefficient", "biot_modulus", "isotropic_permeability"};
        _data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 2), // solid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // biot_coefficient
            pylith::topology::Field::Discretization(0, 2), // biot_modulus
            pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        CPPUNIT_ASSERT(!_data->rheology);
        _data->rheology = new pylith::materials::IsotropicLinearPoroelasticity();CPPUNIT_ASSERT(_data->rheology);
        CPPUNIT_ASSERT(_data->matAuxDB);
        _data->matAuxDB->addValue("solid_density", solid_density, solid_density_units());
        _data->matAuxDB->addValue("fluid_density", fluid_density, fluid_density_units());
        _data->matAuxDB->addValue("fluid_viscosity", fluid_viscosity, fluid_viscosity_units());
        _data->matAuxDB->addValue("porosity", porosity, porosity_units());
        _data->matAuxDB->addValue("shear_modulus", shear_modulus, shear_modulus_units());
        _data->matAuxDB->addValue("drained_bulk_modulus", drained_bulk_modulus, drained_bulk_modulus_units());
        _data->matAuxDB->addValue("biot_coefficient", biot_coefficient, biot_coefficient_units());
        _data->matAuxDB->addValue("fluid_bulk_modulus", fluid_bulk_modulus, fluid_bulk_modulus_units());
        _data->matAuxDB->addValue("solid_bulk_modulus", solid_bulk_modulus, solid_bulk_modulus_units());
        _data->matAuxDB->addValue("isotropic_permeability", isotropic_permeability, isotropic_permeability_units());
        _data->matAuxDB->setCoordSys(*_data->cs);

        _data->faultNumAuxSubfields = 7;
        static const char* _faultAuxSubfields[7] = { "thickness", "porosity", "beta_p", "beta_sigma", "fault_permeability", "fluid_viscosity",  "slip" };
        _data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
            pylith::topology::Field::Discretization(0, 2), // thickness
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // beta_p
            pylith::topology::Field::Discretization(0, 2), // beta_sigma
            pylith::topology::Field::Discretization(0, 2), // fault_permeability
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations =
            const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Fault DB
        CPPUNIT_ASSERT(_data->faultAuxDB);
        _data->faultAuxDB->addValue("thickness", thickness, thickness_units());
        _data->faultAuxDB->addValue("porosity", fault_porosity, fault_porosity_units());
        _data->faultAuxDB->addValue("beta_p", beta_p, beta_p_units());
        _data->faultAuxDB->addValue("beta_sigma", beta_sigma, beta_sigma_units());
        _data->faultAuxDB->addValue("fault_permeability_xx", fault_permeability_xx, fault_permeability_units());
        _data->faultAuxDB->addValue("fault_permeability_yy", fault_permeability_yy, fault_permeability_units());
        _data->faultAuxDB->addValue("fault_permeability_zz", fault_permeability_zz, fault_permeability_units());
        _data->faultAuxDB->addValue("fault_permeability_xy", fault_permeability_xy, fault_permeability_units());
        _data->faultAuxDB->addValue("fluid_viscosity", fault_fluid_viscosity, fault_fluid_viscosity_units());
        _data->faultAuxDB->setCoordSys(*_data->cs);

        // Rupture DB
        CPPUNIT_ASSERT(!_data->kinsrcporo);
        _data->kinsrcporo = new pylith::faults::KinSrcPoroStep();CPPUNIT_ASSERT(_data->kinsrcporo);
        _data->kinsrcporo->originTime(0.0);
        CPPUNIT_ASSERT(_data->ruptureAuxDB);
        _data->ruptureAuxDB->addValue("initiation_time", initiation_time, time_units());
        _data->ruptureAuxDB->addValue("final_slip_opening", finalslip_opening, slip_units());
        _data->ruptureAuxDB->addValue("final_slip_left_lateral", finalslip_leftlateral, slip_units());
        _data->ruptureAuxDB->setCoordSys(*_data->cs);

        // Single aux field for testing
        // _data->faultNumAuxSubfields = 1;
        // static const char* _faultAuxSubfields[1] = { "slip" };
        // _data->faultAuxSubfields = _faultAuxSubfields;
        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
        //     pylith::topology::Field::Discretization(0, 1), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Body MMS Kernels
        std::vector<ResidualKernels> mms_forcing_kernels;
        mms_forcing_kernels.resize(3);
        mms_forcing_kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, NULL);
        mms_forcing_kernels[1] = ResidualKernels("pressure", pylith::feassemble::Integrator::LHS, f0p, NULL);
        mms_forcing_kernels[2] = ResidualKernels("trace_strain", pylith::feassemble::Integrator::LHS, NULL, NULL);

        // Materials
        // _materials.resize(3);
        // { // xneg
        //     pylith::materials::Poroelasticity* material = new pylith::materials::Poroelasticity();assert(material);
        //     material->setFormulation(pylith::problems::Physics::QUASISTATIC);
        //     material->useBodyForce(false);
        //     material->setDescription("Isotropic Linear Poroelasticity Plane Strain");
        //     material->setLabelValue(10);
        //     material->setBulkRheology(_data->rheology);
        //     material->setMMSBodyForceKernels(mms_forcing_kernels);
        //     _materials[0] = material;
        // } // xneg
        _materials.resize(2);
        { // mid
            pylith::materials::Poroelasticity* material = new pylith::materials::Poroelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Poroelasticity Plane Strain");
            material->setLabelValue(20);
            material->setBulkRheology(_data->rheology);
            material->setMMSBodyForceKernels(mms_forcing_kernels);
            _materials[0] = material;
        } // mid
        { // xpos
            pylith::materials::Poroelasticity* material = new pylith::materials::Poroelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescription("Isotropic Linear Poroelasticity Plane Strain");
            material->setLabelValue(15);
            material->setBulkRheology(_data->rheology);
            material->setMMSBodyForceKernels(mms_forcing_kernels);
            _materials[1] = material;
        } // xpos

        // Boundary conditions
        static const PylithInt constrainedDOF_disp[2] = {0, 1};
        static const PylithInt numConstrained_disp = 2;
        static const PylithInt constrainedDOF_pres[1] = {0};
        static const PylithInt numConstrained_pres = 1;
        _bcs.resize(8);
        { // boundary_xneg_disp
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_disp, numConstrained_disp);
            bc->setUserFn(solnkernel_disp);
            _bcs[0] = bc;
        } // boundary_xneg_disp
        { // boundary_xpos_disp
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_disp, numConstrained_disp);
            bc->setUserFn(solnkernel_disp);
            _bcs[1] = bc;
        } // boundary_xpos_disp
        { // boundary_yneg_disp
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_disp, numConstrained_disp);
            bc->setUserFn(solnkernel_disp);
            _bcs[2] = bc;
        } // boundary_yneg_disp
        { // boundary_ypos_disp
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("displacement");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_disp, numConstrained_disp);
            bc->setUserFn(solnkernel_disp);
            _bcs[3] = bc;
        } // boundary_ypos_disp
        { // boundary_xneg_pres
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_pres, numConstrained_pres);
            bc->setUserFn(solnkernel_pres);
            _bcs[4] = bc;
        } // boundary_xneg_pres
        { // boundary_xpos_pres
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_xpos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_pres, numConstrained_pres);
            bc->setUserFn(solnkernel_pres);
            _bcs[5] = bc;
        } // boundary_xpos_disp
        { // boundary_yneg_disp
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_yneg");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_pres, numConstrained_pres);
            bc->setUserFn(solnkernel_pres);
            _bcs[6] = bc;
        } // boundary_yneg_pres
        { // boundary_ypos_pres
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setSubfieldName("pressure");
            bc->setLabelName("boundary_ypos");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF_pres, numConstrained_pres);
            bc->setUserFn(solnkernel_pres);
            _bcs[7] = bc;
        } // boundary_ypos_pres

        // Faults
        _faults.resize(1);
        { // xpos
            pylith::faults::FaultCohesiveKinPoro* fault = new pylith::faults::FaultCohesiveKinPoro();
            fault->setCohesiveLabelValue(100);
            fault->setSurfaceLabelName("fault_xpos");

            const int numRuptures = 1;
            const char* ruptureNames[1] = { "rupture" };
            pylith::faults::KinSrcPoro* ruptures[1] = { _data->kinsrcporo };
            fault->setEqRuptures(ruptureNames, numRuptures, ruptures, numRuptures);
            _faults[0] = fault;
        } // xpos
        pylith::utils::PetscOptions options;
        options.add("-pc_type", "jacobi");
        options.add("-fieldsplit_displacement_pc_type", "lu");
        options.override ();
    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        const pylith::topology::Field* solution = _problem->getSolution();
        CPPUNIT_ASSERT(solution);

        PetscDM dm = solution->getDM();
        PetscDMLabel label;
        PetscIS is;
        PetscInt cohesiveCell;
        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(dm, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pres, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_tracestrain, dm);CPPUNIT_ASSERT(!err);
        err = DMGetLabel(dm, pylith::topology::Mesh::cells_label_name, &label);CPPUNIT_ASSERT(!err);
        err = DMLabelGetStratumIS(label, _faults[0]->getCohesiveLabelValue(), &is);CPPUNIT_ASSERT(!err);
        err = ISGetMinMax(is, &cohesiveCell, NULL);CPPUNIT_ASSERT(!err);
        err = ISDestroy(&is);CPPUNIT_ASSERT(!err);
        err = DMGetCellDS(dm, cohesiveCell, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pres, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_tracestrain, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 3, solnkernel_lagrangemultiplier, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 4, solnkernel_faultpressure, dm);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestFaultKinPoro2D_OneFaultSlip

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP1 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_TriP1,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);
        _allowZeroResidual = true;

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // pressure
            pylith::topology::Field::Discretization(1, 2), // trace_strain
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(1, 2, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP2 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_TriP2,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 2), // solid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // biot_coefficient
            pylith::topology::Field::Discretization(0, 2), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 2), // thickness
        //     pylith::topology::Field::Discretization(0, 2), // porosity
        //     pylith::topology::Field::Discretization(0, 2), // beta_p
        //     pylith::topology::Field::Discretization(0, 2), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 2), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 2), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(3, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // pressure
            pylith::topology::Field::Discretization(2, 2), // trace_strain
            pylith::topology::Field::Discretization(3, 2, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP3 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_TriP3,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 3), // solid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 3), // porosity
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // biot_coefficient
            pylith::topology::Field::Discretization(0, 3), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 3), // thickness
        //     pylith::topology::Field::Discretization(0, 3), // porosity
        //     pylith::topology::Field::Discretization(0, 3), // beta_p
        //     pylith::topology::Field::Discretization(0, 3), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 3), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 3), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(4, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // pressure
            pylith::topology::Field::Discretization(3, 3), // trace_strain
            pylith::topology::Field::Discretization(4, 3, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP4 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_TriP4,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 4), // solid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // biot_coefficient
            pylith::topology::Field::Discretization(0, 4), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 4), // thickness
        //     pylith::topology::Field::Discretization(0, 4), // porosity
        //     pylith::topology::Field::Discretization(0, 4), // beta_p
        //     pylith::topology::Field::Discretization(0, 4), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 4), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 4), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(5, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // pressure
            pylith::topology::Field::Discretization(4, 4), // trace_strain
            pylith::topology::Field::Discretization(5, 4, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ1 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_QuadQ1,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 1), // disp
            pylith::topology::Field::Discretization(1, 1), // pressure
            pylith::topology::Field::Discretization(1, 1), // trace_strain
            pylith::topology::Field::Discretization(2, 1, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(1, 1, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ2 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_QuadQ2,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 2), // solid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // biot_coefficient
            pylith::topology::Field::Discretization(0, 2), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 2), // thickness
        //     pylith::topology::Field::Discretization(0, 2), // porosity
        //     pylith::topology::Field::Discretization(0, 2), // beta_p
        //     pylith::topology::Field::Discretization(0, 2), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 2), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 2), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(3, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // pressure
            pylith::topology::Field::Discretization(2, 2), // trace_strain
            pylith::topology::Field::Discretization(3, 2, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ3 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_QuadQ3,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 3), // solid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 3), // porosity
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // biot_coefficient
            pylith::topology::Field::Discretization(0, 3), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 3), // thickness
        //     pylith::topology::Field::Discretization(0, 3), // porosity
        //     pylith::topology::Field::Discretization(0, 3), // beta_p
        //     pylith::topology::Field::Discretization(0, 3), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 3), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 3), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(4, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // pressure
            pylith::topology::Field::Discretization(3, 3), // trace_strain
            pylith::topology::Field::Discretization(4, 3, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ4 :
    public pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_OneFaultSlip_QuadQ4,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_OneFaultSlip::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 4), // solid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_density
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // biot_coefficient
            pylith::topology::Field::Discretization(0, 4), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 4), // isotropic_permeability
        };
        _data->matAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_matAuxDiscretizations);

        // static const pylith::topology::Field::Discretization _faultAuxDiscretizations[7] = {
        //     pylith::topology::Field::Discretization(0, 4), // thickness
        //     pylith::topology::Field::Discretization(0, 4), // porosity
        //     pylith::topology::Field::Discretization(0, 4), // beta_p
        //     pylith::topology::Field::Discretization(0, 4), // beta_sigma
        //     pylith::topology::Field::Discretization(0, 4), // fault_permeability
        //     pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
        //     pylith::topology::Field::Discretization(0, 4), // slip
        // };
        // _data->faultAuxDiscretizations =
        // const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(5, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // pressure
            pylith::topology::Field::Discretization(4, 4), // trace_strain
            pylith::topology::Field::Discretization(5, 4, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_OneFaultSlip_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_OneFaultSlip_QuadQ4);

// End of file
