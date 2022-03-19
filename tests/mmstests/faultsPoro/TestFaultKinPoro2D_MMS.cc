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

/** @file tests/mmstests/faults/TestFaultKinPoro2D_MMS.cc
 *
 * Square domain of sides 8.0 km with a through-going fault running
 * through the center in the y-direction. The two opposing sides each
 * move as rigid blocks with 3.0 m of right-lateral slip.
 */

#include <portinfo>

#include "TestFaultKinPoro.hh" // ISA TestFaultKinPoro

#include "pylith/faults/FaultCohesiveKinPoro.hh" // USES FaultCohesiveKinPoro
#include "pylith/faults/KinSrcPoroStep.hh" // USES KinSrcPoroStep
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestFaultKinPoro2D_MMS;

        class TestFaultKinPoro2D_MMS_TriP1;
        class TestFaultKinPoro2D_MMS_TriP2;
        class TestFaultKinPoro2D_MMS_TriP3;
        class TestFaultKinPoro2D_MMS_TriP4;

        class TestFaultKinPoro2D_MMS_QuadQ1;
        class TestFaultKinPoro2D_MMS_QuadQ2;
        class TestFaultKinPoro2D_MMS_QuadQ3;
        class TestFaultKinPoro2D_MMS_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS :
    public pylith::mmstests::TestFaultKinPoro {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

    // Solid Density
    static double solid_density(const double x,
                                const double y) {
        return 2500.0;
    } // solid_density

    static const char* solid_density_units(void) {
        return "kg/m**3";
    } // solid_density_units

    // Fluid Density
    static double fluid_density(const double x,
                                const double y) {
        return 1000.0;
    } // fluid_density

    static const char* fluid_density_units(void) {
        return "kg/m**3";
    } // fluid_density_units

    // Fluid Viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 1.0;
    } // fluid_viscosity

    static const char* fluid_viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.5;
    } // porosity

    static const char* porosity_units(void) {
        return "none";
    } // porosity_units

    // Shear Modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 0.5;
    } // shear_modulus

    static const char* shear_modulus_units(void) {
        return "Pa";
    } // shear_modulus_units

    // Drained Bulk Modulus
    static double drained_bulk_modulus(const double x,
                                       const double y) {
        return 0.5;
    } // shear_modulus

    static const char* drained_bulk_modulus_units(void) {
        return "Pa";
    } // drained_bulk_modulus_units

    // Biot Coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 1.0;
    } // alpha

    static const char* biot_coefficient_units(void) {
        return "none";
    } // alpha_units

    // Fluid Bulk Modulus
    static double fluid_bulk_modulus(const double x,
                                     const double y) {
        return 0.5;
    } // fluid_bulk_modulus

    static const char* fluid_bulk_modulus_units(void) {
        return "Pa";
    } // fluid_bulk_modulus_units

    // Solid Bulk Modulus
    static double solid_bulk_modulus(const double x,
                                     const double y) {
        return 0.5;
    } // solid_bulk_modulus

    static const char* solid_bulk_modulus_units(void) {
        return "Pa";
    } // solid_bulk_modulus_units

    // Isotropic Permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 0.5;
    } // isotropic_permeability

    static const char* isotropic_permeability_units(void) {
        return "m**2";
    } // isotropic_permeability_units

    // Poroelastic fault auxiliary components

    // Thickness
    static double thickness(const double x,
                            const double y) {
        return 1.0;
    } // thickness

    static const char* thickness_units(void) {
        return "m";
    } // thickness_units

    // Fault porosity
    static double fault_porosity(const double x,
                                 const double y) {
        return 0.5;
    } // fault_porosity

    static const char* fault_porosity_units(void) {
        return "one";
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

    // Tangential Permeability
    static double tangential_permeability(const double x,
                                          const double y) {
        return 1.0;
    } // tangential_permeability

    static const char* tangential_permeability_units(void) {
        return "m**2";
    } // tangential_permeability_units

    // Normal Permeability
    static double normal_permeability(const double x,
                                      const double y) {
        return 1.0;
    } // normal_permeability

    static const char* normal_permeability_units(void) {
        return "m**2";
    } // normal_permeability_units

    // Fault Fluid Viscosity
    static double fault_fluid_viscosity(const double x,
                                        const double y) {
        return 1.0;
    } // fault_fluid_viscosity

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
        return -1.5;
    } // slip_leftlateral

    static const char* slip_units(void) {
        return "m";
    } // slip_units

    // Solution subfields.

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        if (y > 0.0) {
            return 1.0 - y*y;
        } else if (y < 0.0) {
            return y*y - 1.0;
        } else {
            return 0.0;
        }
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        return y*abs(y);
    } // disp_y

    // Volumetric Strain
    static double trace_strain(const double x,
                               const double y,
                               const double t) {
        if (y > 0) {
            return -2.0*y;
        } else if (y < 0) {
            return 2.0*y;
        } else {
            return 0.0;
        } // if/else
    } // trace_strain

    // Pressure
    static double pressure(const double x,
                           const double y,
                           const double t) {
        return t*(y*y - 1);
    }

    // Fault Traction, Normal
    static double faulttraction_n(const double x,
                                  const double y,
                                  const double t) {
        return 0.0;
    } // faulttraction_n

    // Fault Traction, Tangential
    static double faulttraction_t(const double x,
                                  const double y,
                                  const double t) {
        return t;
    } // faulttraction_t

    // Fault Pressure
    static double faultpressure(const double x,
                                const double y,
                                const double t) {
        return -4.0*y + t;
    } // faultpressure

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

        s[0] = disp_x(x[0], x[1], t);
        s[1] = disp_y(x[0], x[1], t);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_pressure(PetscInt spaceDim,
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
    } // solnkernel_pressure

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

        s[0] = faulttraction_n(x[0], x[1], t);
        s[1] = faulttraction_t(x[0], x[1], t);

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

        s[0] = faultpressure(x[0], x[1], t);

        return 0;
    } // solnkernel_faultpressure

protected:

    void setUp(void) {
        TestFaultKinPoro::setUp();

        // Overwrite component name for control of journals at test level.
        GenericComponent::setName("TestFaultKinPoro2D_MMS");

        CPPUNIT_ASSERT(!_data);
        _data = new TestFaultKinPoro_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->setLengthScale(1.0e+03);
        _data->normalizer->setTimeScale(2.0);
        _data->normalizer->setPressureScale(2.25e+10);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->matNumAuxSubfields = 10;
        static const char* _matAuxSubfields[10] = {"solid_density", "fluid_density", "fluid_viscosity", "porosity", "shear_modulus", "drained_bulk_modulus", "biot_coefficient", "fluid_bulk_modulus", "solid_bulk_modulus", "isotropic_permeability"};
        _data->matAuxSubfields = _matAuxSubfields;
        static const pylith::topology::Field::Discretization _matAuxDiscretizations[10] = {
            pylith::topology::Field::Discretization(0, 1), // solid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // drained_bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // biot_coefficient
            pylith::topology::Field::Discretization(0, 1), // fluid_bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // solid_bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // isotropic_permeability
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

        CPPUNIT_ASSERT(!_data->kinsrcporo);
        _data->kinsrcporo = new pylith::faults::KinSrcPoroStep();CPPUNIT_ASSERT(_data->kinsrcporo);
        _data->kinsrcporo->originTime(0.0);
        CPPUNIT_ASSERT(_data->faultAuxDB);
        _data->faultAuxDB->addValue("thickness", thickness, thickness_units());
        _data->faultAuxDB->addValue("porosity", porosity, porosity_units());
        _data->faultAuxDB->addValue("beta_p", beta_p, beta_p_units());
        _data->faultAuxDB->addValue("beta_sigma", beta_sigma, beta_sigma_units());
        _data->faultAuxDB->addValue("permeability_normal", normal_permeability, normal_permeability_units());
        _data->faultAuxDB->addValue("permeability_tangential", tangential_permeability, tangential_permeability_units());
        _data->faultAuxDB->addValue("fluid_viscosity", fluid_viscosity, fluid_viscosity_units());
        _data->faultAuxDB->addValue("initiation_time", initiation_time, time_units());
        _data->faultAuxDB->addValue("final_slip_opening", finalslip_opening, slip_units());
        _data->faultAuxDB->addValue("final_slip_left_lateral", finalslip_leftlateral, slip_units());
        _data->faultAuxDB->setCoordSys(*_data->cs);

        _data->faultNumAuxSubfields = 8;
        static const char* _faultAuxSubfields[8] = {"thickness", "porosity", "beta_p", "beta_sigma", "permeability_tangential", "permeability_normal", "fluid_viscosity", "slip" };
        _data->faultAuxSubfields = _faultAuxSubfields;
        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 1), // thickness
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // beta_p
            pylith::topology::Field::Discretization(0, 1), // beta_sigma
            pylith::topology::Field::Discretization(0, 1), // permeability_tangential
            pylith::topology::Field::Discretization(0, 1), // permeability_normal
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        // Materials
        _materials.resize(2);
        { // xneg
            pylith::materials::Poroelasticity* material = new pylith::materials::Poroelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescriptiveLabel("Isotropic Linear Poroelasticity Plane Strain");
            material->setMaterialId(10);
            material->setBulkRheology(_data->rheology);
            _materials[0] = material;
        } // xneg
        { // xpos
            pylith::materials::Poroelasticity* material = new pylith::materials::Poroelasticity();assert(material);
            material->setFormulation(pylith::problems::Physics::QUASISTATIC);
            material->useBodyForce(false);
            material->setDescriptiveLabel("Isotropic Linear Poroelasticity Plane Strain");
            material->setMaterialId(20);
            material->setBulkRheology(_data->rheology);
            _materials[1] = material;
        } // xpos

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        static const PylithInt constrainedDOF_pressure[1] = {0};
        static const PylithInt numConstrained_pressure = 1;
        _bcs.resize(8);
        { // boundary_xpos
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_xpos");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[0] = bc;
        } // boundary_xpos
        { // boundary_xneg
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_xneg");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[1] = bc;
        } // boundary_xneg
        { // boundary_xpos_neu
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF_pressure, numConstrained_pressure);
            bc->setMarkerLabel("boundary_xpos_neu");
            bc->setSubfieldName("pressure");
            bc->setUserFn(solnkernel_pressure);
            _bcs[2] = bc;
        } // boundary_xpos_neu
        { // boundary_xneg_neu
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF_pressure, numConstrained_pressure);
            bc->setMarkerLabel("boundary_xneg_neu");
            bc->setSubfieldName("pressure");
            bc->setUserFn(solnkernel_pressure);
            _bcs[3] = bc;
        } // boundary_xneg_neu
        { // boundary_ypos
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_ypos");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[4] = bc;
        } // boundary_ypos
        { // boundary_yneg
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF, numConstrained);
            bc->setMarkerLabel("boundary_yneg");
            bc->setSubfieldName("displacement");
            bc->setUserFn(solnkernel_disp);
            _bcs[5] = bc;
        } // boundary_yneg
        { // boundary_ypos_neu
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF_pressure, numConstrained_pressure);
            bc->setMarkerLabel("boundary_ypos_neu");
            bc->setSubfieldName("pressure");
            bc->setUserFn(solnkernel_pressure);
            _bcs[6] = bc;
        } // boundary_ypos_neu
        { // boundary_yneg_neu
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();
            bc->setConstrainedDOF(constrainedDOF_pressure, numConstrained_pressure);
            bc->setMarkerLabel("boundary_yneg_neu");
            bc->setSubfieldName("pressure");
            bc->setUserFn(solnkernel_pressure);
            _bcs[7] = bc;
        } // boundary_yneg_neu

        _fault->setInterfaceId(100);
        _fault->setSurfaceMarkerLabel("fault");

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscDM dm = _solution->getDM();
        PetscDMLabel label;
        PetscIS is;
        PetscInt cohesiveCell;
        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(dm, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, dm);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_tracestrain, dm);CPPUNIT_ASSERT(!err);
        err = DMGetLabel(dm, "material-id", &label);CPPUNIT_ASSERT(!err);
        err = DMLabelGetStratumIS(label, _fault->getInterfaceId(), &is);CPPUNIT_ASSERT(!err);
        err = ISGetMinMax(is, &cohesiveCell, NULL);CPPUNIT_ASSERT(!err);
        err = ISDestroy(&is);CPPUNIT_ASSERT(!err);
        err = DMGetCellDS(dm, cohesiveCell, &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_tracestrain, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 3, solnkernel_lagrangemultiplier, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 4, solnkernel_faultpressure, NULL);CPPUNIT_ASSERT(!err);

    } // _setExactSolution

}; // TestFaultKinPoro2D_MMS

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_TriP1 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_TriP1,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
        CPPUNIT_ASSERT(_data);
        _allowZeroResidual = true;

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 1), // disp
            pylith::topology::Field::Discretization(2, 1), // pressure
            pylith::topology::Field::Discretization(2, 1), // trace_strain
            pylith::topology::Field::Discretization(2, 1, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 1, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_TriP2 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_TriP2,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 2), // thickness
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // beta_p
            pylith::topology::Field::Discretization(0, 2), // beta_sigma
            pylith::topology::Field::Discretization(0, 2), // permeability_tangential
            pylith::topology::Field::Discretization(0, 2), // permeability_normal
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // pressure
            pylith::topology::Field::Discretization(2, 2), // trace_strain
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_TriP3 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_TriP3,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 3), // thickness
            pylith::topology::Field::Discretization(0, 3), // porosity
            pylith::topology::Field::Discretization(0, 3), // beta_p
            pylith::topology::Field::Discretization(0, 3), // beta_sigma
            pylith::topology::Field::Discretization(0, 3), // permeability_tangential
            pylith::topology::Field::Discretization(0, 3), // permeability_normal
            pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // pressure
            pylith::topology::Field::Discretization(3, 3), // trace_strain
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_TriP4 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_TriP4,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 4), // thickness
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // beta_p
            pylith::topology::Field::Discretization(0, 4), // beta_sigma
            pylith::topology::Field::Discretization(0, 4), // permeability_tangential
            pylith::topology::Field::Discretization(0, 4), // permeability_normal
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // pressure
            pylith::topology::Field::Discretization(4, 4), // trace_strain
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_TriP4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ1 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_QuadQ1,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 1), // disp
            pylith::topology::Field::Discretization(2, 1), // pressure
            pylith::topology::Field::Discretization(2, 1), // trace_strain
            pylith::topology::Field::Discretization(2, 1, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 1, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ2 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_QuadQ2,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 2), // thickness
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // beta_p
            pylith::topology::Field::Discretization(0, 2), // beta_sigma
            pylith::topology::Field::Discretization(0, 2), // permeability_tangential
            pylith::topology::Field::Discretization(0, 2), // permeability_normal
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(2, 2), // pressure
            pylith::topology::Field::Discretization(2, 2), // trace_strain
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(2, 2, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ3 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_QuadQ3,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 3), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(3, 3), // pressure
            pylith::topology::Field::Discretization(3, 3), // trace_strain
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(3, 3, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ4 :
    public pylith::mmstests::TestFaultKinPoro2D_MMS {
    CPPUNIT_TEST_SUB_SUITE(TestFaultKinPoro2D_MMS_QuadQ4,
                           TestFaultKinPoro);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestFaultKinPoro2D_MMS::setUp();
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

        static const pylith::topology::Field::Discretization _faultAuxDiscretizations[8] = {
            pylith::topology::Field::Discretization(0, 4), // thickness
            pylith::topology::Field::Discretization(0, 4), // porosity
            pylith::topology::Field::Discretization(0, 4), // beta_p
            pylith::topology::Field::Discretization(0, 4), // beta_sigma
            pylith::topology::Field::Discretization(0, 4), // permeability_tangential
            pylith::topology::Field::Discretization(0, 4), // permeability_normal
            pylith::topology::Field::Discretization(0, 4), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 4), // slip
        };
        _data->faultAuxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_faultAuxDiscretizations);

        _data->numSolnSubfields = 5;
        static const pylith::topology::Field::Discretization _solnDiscretizations[5] = {
            pylith::topology::Field::Discretization(4, 4), // disp
            pylith::topology::Field::Discretization(4, 4), // pressure
            pylith::topology::Field::Discretization(4, 4), // trace_strain
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // lagrange_multiplier_fault
            pylith::topology::Field::Discretization(4, 4, 1, -1, true), // fault_pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // TestFaultKinPoro2D_MMS_QuadQ4
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestFaultKinPoro2D_MMS_QuadQ4);

// End of file
