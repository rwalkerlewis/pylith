// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestIsotropicLinearPoroelasticity.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn
#include "pylith/bc/NeumannTimeDependent.hh" // USES NeumannTimeDependent

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#define PI 3.14159265


namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearPoroelasticity_Terzaghi;

        class TestIsotropicLinearPoroelasticity_Terzaghi_TriP1;
        class TestIsotropicLinearPoroelasticity_Terzaghi_TriP2;
        class TestIsotropicLinearPoroelasticity_Terzaghi_TriP3;
        //class TestIsotropicLinearPoroelasticity_Terzaghi_TriP4;

        class TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ2;
        class TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ3;
        //class TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity {
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double GACC;
    static const double YMAX;
    static const double YMIN;
    static const double EXTERNALFORCE;
    static const int NITERATIONS;

    /// Spatial database user functions for auxiliary subfields (includes derived fields).

    // Porosity
    static double porosity(const double x,
                          const double y) {
        return 0.19;
    } // porosity

    static const char* porosity_units(void) {
        return "dimensionless";
    } // porosity_units

    // Solid Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Fluid density
    static double fluid_density(const double x,
                                const double y) {
        return 1000.0;
    } // fluid_density

    static const char* fluid_density_units(void) {
        return "kg/m**3";
    } // fluid_density_units

    // Fluid viscosity
    static double fluid_viscosity(const double x,
                                  const double y) {
        return 0.001;
    } // fluid_viscosity

    static const char* fluid_viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    // Shear modulus
    static double shear_modulus(const double x,
                                const double y) {
        return 6e9;
    } // shear_modulus

    static const char* shear_modulus_units(void) {
        return "Pa";
    } // shear_modulus_units

    // Solid Bulk Modulus
    static double bulk_modulus(const double x,
                               const double y) {
        return 3.6e10;
    } // bulk_modulus

    static const char* bulk_modulus_units(void) {
        return "Pa";
    } // bulk_modulus_units

    // Biot coefficient
    static double biot_coefficient(const double x,
                                   const double y) {
        return 0.778;
    } // biot_coefficient

    static const char* biot_coefficient_units(void) {
        return "dimensionless";
    } // biot_coefficient_units

    // Isotropic permeability
    static double isotropic_permeability(const double x,
                                         const double y) {
        return 1.90e-13;
    } // isotropic_permeability

    static const char* isotropic_permeability_units(void) {
        return "m**2";
    } // isotropic_permeability_units

    // Fluid Bulk Modulus
    static double fluid_bulk_modulus(const double x,
                               const double y) {
        return 2e9;
    } // fluid_bulk_modulus

    static const char* fluid_bulk_modulus_units(void) {
        return "Pa";
    } // fluid_bulk_modulus_units

    static double gravityAcc_x(const double x,
                               const double y) {
        return 0.0;
    } // gravityAcc_x

    static double gravityAcc_y(const double x,
                               const double y) {
        return 0.0;
    } // gravityAcc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

    // Derived fields

    // Biot Modulus
    static double biot_modulus(const double x,
                              const double y){
        return bulk_modulus(x,y) / (biot_coefficient(x,y) - porosity(x,y)*(1.0 - bulk_modulus(x,y)/fluid_bulk_modulus(x,y)));
    } // biot_modulus

    static const char* biot_modulus_units(void) {
        return "Pa";
    } // biot_modulus_units

    // Drained Bulk Modulus
    static double drained_bulk_modulus(const double x,
                                      const double y){
        return bulk_modulus(x,y) * (1.0 - biot_coefficient(x,y));
    } // drained_bulk_modulus

    static const char* drained_bulk_modulus_units(void) {
        return "Pa";
    } // drained_bulk_modulus_units

    // Undrained Bulk Modulus
    static double undrained_bulk_modulus(const double x,
                                        const double y){
        return drained_bulk_modulus(x,y) + biot_coefficient(x,y)*biot_coefficient(x,y)*biot_modulus(x,y);
    } // undrained_bulk_modulus

    static const char* undrained_bulk_modulus_units(void) {
        return "Pa";
    } // undrained_bulk_modulus_units

    // Poisson's Ratio
    static double poisson_ratio(const double x,
                               const double y){
        return (3.0*drained_bulk_modulus(x,y) - 2.0*shear_modulus(x,y)) / (2.0 * (2.0*drained_bulk_modulus(x,y) + shear_modulus(x,y)));
    } // poisson_ratio

    static const char* poisson_ratio_units(void) {
        return "dimensionless";
    } // poisson_ratio_units

    // Undrained Poisson's Ratio
    static double undrained_poisson_ratio(const double x,
                               const double y){
        return (3.0*undrained_bulk_modulus(x,y) - 2.0*shear_modulus(x,y)) / (2.0 * (2.0*undrained_bulk_modulus(x,y) + shear_modulus(x,y)));
    } // undrained_poisson_ratio

    static const char* undrained_poisson_ratio_units(void) {
        return "dimensionless";
    } // undrained_poisson_ratio_units

    // Consolidation Coefficient
    static double consolidation_coefficient(const double x,
                                           const double y){
        return (isotropic_permeability(x,y) / fluid_viscosity(x,y)) / (porosity(x,y)/fluid_bulk_modulus(x,y) + 3.0/(3.0*drained_bulk_modulus(x,y) + 2.0*shear_modulus(x,y)));
    } // consolidation_coefficient

    static const char* consolidation_coefficient_units(void) {
        return "m**2 / s";
    } // consolidation_coefficient_units

    // Neumann Boundary Condition

    // Initial Amplitude Tangential
    static double initial_amplitude_tangential(const double x,
                                           const double y){
        return 0.0;
    } // initial_amplitude_tangential

    static const char* initial_amplitude_tangential_units(void) {
        return "Pa";
    } // initial_amplitude_tangential_units

    // Initial Amplitude Normal
    static double initial_amplitude_normal(const double x,
                                           const double y){
        return -1000.0;
    } // initial_amplitude_normal

    static const char* initial_amplitude_normal_units(void) {
        return "Pa";
    } // initial_amplitude_normal_units

    // Solution fields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {


       return 0.0 / LENGTHSCALE;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {

        const double L = YMAX - YMIN;
        const double zstar = y /  L;
        const double tstar = (consolidation_coefficient(x,y) * t) / (4.0 * L*L);

        // Series term
        double F2 = 0.0;
        for (int m = 1; m < NITERATIONS*2+1; m++) {
          if (m%2 == 1) {
            F2 += 8.0/(m*m * PI*PI) * cos( (m*PI*zstar) / 2.0 ) * (1.0 - exp(-1*m*m*PI*PI*tstar) );
          }
        }
        double disp = ( EXTERNALFORCE * L * (1.0 - 2.0*undrained_poisson_ratio(x,y)) )/( 2.0*shear_modulus(x,y)*(1.0 - undrained_poisson_ratio(x,y)) ) * (1.0 - zstar) + ( EXTERNALFORCE * L * (undrained_poisson_ratio(x,y) - poisson_ratio(x,y)) )/( 2.0*shear_modulus(x,y)*(1.0 - undrained_poisson_ratio(x,y))*(1.0 - poisson_ratio(x,y)))  * F2;
        return disp / LENGTHSCALE;
    } // disp_y

    // Pressure
    static double pressure(const double x,
                           const double y,
                           const double t) {

        const double L = YMAX - YMIN;
        const double zstar = y /  L;
        const double tstar = (consolidation_coefficient(x,y) * t) / (4.0 * L*L);

        const double eta = biot_coefficient(x,y) * ( 1.0 - 2.0*poisson_ratio(x,y)) / (2.0*(1.0 - poisson_ratio(x,y)));
        const double S = porosity(x,y)/fluid_bulk_modulus(x,y) + 3.0/(3.0*drained_bulk_modulus(x,y) + 2.0*shear_modulus(x,y));

        // Series term
        double F1 = 0.0;
        for (int m = 1; m < NITERATIONS*2+1; m++) {
          if (m%2 == 1) {
            F1 += 4.0/(m * PI) * sin( (m*PI*zstar) / 2.0 ) * exp(-1*m*m*PI*PI*tstar);
          }
        }
        double pressure = (EXTERNALFORCE*eta)/(shear_modulus(x,y)*S) * F1;

        return pressure / PRESSURESCALE;
    } // pressure

    // Trace Strain
    static double trace_strain(const double x,
                               const double y,
                               const double t) {

        const double L = YMAX - YMIN;
        const double zstar = y /  L;
        const double tstar = (consolidation_coefficient(x,y) * t) / (4.0 * L*L);

        double C = ( EXTERNALFORCE * L * (1.0 - 2.0*undrained_poisson_ratio(x,y)) )/( 2.0*shear_modulus(x,y)*(1.0 - undrained_poisson_ratio(x,y)) );
        double D = ( EXTERNALFORCE * L * (undrained_poisson_ratio(x,y) - poisson_ratio(x,y)) )/( 2.0*shear_modulus(x,y)*(1.0 - undrained_poisson_ratio(x,y))*(1.0 - poisson_ratio(x,y)));

        // Series term
        double dF2_dzstar = 0.0;
        for (int m = 1; m < NITERATIONS*2+1; m++) {
          if (m%2 == 1) {
            dF2_dzstar += -4.0/(m * PI) * sin( (m*PI*zstar) / 2.0 ) * (1.0 - exp(-1*m*m*PI*PI*tstar) );
          }
        }

        return -1.0*C / L + (D/L)*dF2_dzstar;
      } // trace_strain

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);
        CPPUNIT_ASSERT(x);



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

    static PetscErrorCode solnkernel_trace_strain(PetscInt spaceDim,
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
    } // solnkernel_trace_strain

protected:

    void setUp(void) {
        TestIsotropicLinearPoroelasticity::setUp();

        // Overwrite component names for control of debugging info at test level.
        GenericComponent::setName("TestIsotropicLinearPoroelasticity_Terzaghi");
        journal::debug_t debug(GenericComponent::getName());
        //debug.activate(); // DEBUGGING

        CPPUNIT_ASSERT(!_data);
        _data = new TestPoroelasticity_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);
        _data->cs->initialize();

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->lengthScale(LENGTHSCALE);
        _data->normalizer->timeScale(TIMESCALE);
        _data->normalizer->pressureScale(PRESSURESCALE);
        _data->normalizer->computeDensityScale();

        _data->startTime = 0.0;
        _data->endTime = 5.0;
        _data->timeStep = 1.0;

        // solnDiscretizations set in derived class.

        // Material Auxililary
        _data->numAuxSubfields = 9;
        static const char* _auxSubfields[9] = {
            "porosity",
            "density",
            "fluid_density",
            "fluid_viscosity",
            "shear_modulus",
            "bulk_modulus",
            "biot_coefficient",
            "isotropic_permeability",
            "fluid_bulk_modulus",
        };
        _data->auxSubfields = _auxSubfields;

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // biot_coefficient
            pylith::topology::Field::Discretization(0, 1), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 1), // fluid_bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("porosity", porosity, porosity_units());
        _data->auxDB->addValue("density", density, density_units());
        _data->auxDB->addValue("fluid_density", fluid_density, fluid_density_units());
        _data->auxDB->addValue("fluid_viscosity", fluid_viscosity, fluid_viscosity_units());
        _data->auxDB->addValue("shear_modulus", shear_modulus, shear_modulus_units());
        _data->auxDB->addValue("bulk_modulus", bulk_modulus, bulk_modulus_units());
        _data->auxDB->addValue("biot_coefficient", biot_coefficient, biot_coefficient_units());
        _data->auxDB->addValue("isotropic_permeability", isotropic_permeability, isotropic_permeability_units());
        _data->auxDB->addValue("fluid_bulk_modulus", fluid_bulk_modulus, fluid_bulk_modulus_units());
        _data->auxDB->coordsys(*_data->cs);

        // Traction DB
        _data->numTractionSubfields = 2;
        static const char* _tractionSubfields[2] = {
            "initial_amplitude_tangential",
            "initial_amplitude_normal",
        };
        _data->tractionSubfields = _tractionSubfields;

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 1), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 1), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

        CPPUNIT_ASSERT(_data->tractionDB);
        _data->tractionDB->addValue("initial_amplitude_tangential", initial_amplitude_tangential, initial_amplitude_tangential_units());
        _data->tractionDB->addValue("initial_amplitude_normal", initial_amplitude_normal, initial_amplitude_normal_units());
        _data->tractionDB->coordsys(*_data->cs);

        // Material flags
        CPPUNIT_ASSERT(_material);
        _material->useBodyForce(false);
        _rheology->useReferenceState(false);

        _material->setDescriptiveLabel("Isotropic Linear Poroelasticity");
        _material->setMaterialId(24);

        // X Negative Displacement
        static const PylithInt constrainedDispDOF_xneg[1] = {0};
        static const PylithInt numConstrainedDisp_xneg = 1;
        _bcDisplacement_xneg->setConstrainedDOF(constrainedDispDOF_xneg, numConstrainedDisp_xneg);
        _bcDisplacement_xneg->setMarkerLabel("xneg");
        _bcDisplacement_xneg->setSubfieldName("displacement");
        _bcDisplacement_xneg->setUserFn(solnkernel_disp);

        // X Positive Displacement
        static const PylithInt constrainedDispDOF_xpos[1] = {0};
        static const PylithInt numConstrainedDisp_xpos = 1;
        _bcDisplacement_xpos->setConstrainedDOF(constrainedDispDOF_xpos, numConstrainedDisp_xpos);
        _bcDisplacement_xpos->setMarkerLabel("xpos");
        _bcDisplacement_xpos->setSubfieldName("displacement");
        _bcDisplacement_xpos->setUserFn(solnkernel_disp);

        // Y Negative Displacement
        static const PylithInt constrainedDispDOF_yneg[1] = {1};
        static const PylithInt numConstrainedDisp_yneg = 1;
        _bcDisplacement_yneg->setConstrainedDOF(constrainedDispDOF_yneg, numConstrainedDisp_yneg);
        _bcDisplacement_yneg->setMarkerLabel("yneg");
        _bcDisplacement_yneg->setSubfieldName("displacement");
        _bcDisplacement_yneg->setUserFn(solnkernel_disp);

        // Y Positive Pressure
        static const PylithInt constrainedPressureDOF_ypos[1] = {0};
        static const PylithInt numConstrainedPressure_ypos = 1;
        _bcPressure_ypos->setConstrainedDOF(constrainedPressureDOF_ypos, numConstrainedPressure_ypos);
        _bcPressure_ypos->setMarkerLabel("ypos");
        _bcPressure_ypos->setSubfieldName("pressure");
        _bcPressure_ypos->setUserFn(solnkernel_pressure);

        // Y Positive Traction
        _bcTraction_ypos->useInitial(true);
        _bcTraction_ypos->useRate(false);
        _bcTraction_ypos->setScaleName("pressure");
        _bcTraction_ypos->setMarkerLabel("ypos");
        _bcTraction_ypos->setSubfieldName("displacement");

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 2, solnkernel_trace_strain, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestIsotropicLinearPoroelasticity_Terzaghi
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::LENGTHSCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::TIMESCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::PRESSURESCALE = 1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::GACC = 9.80665;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::YMAX = +1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::YMIN = -1.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::EXTERNALFORCE = 1000;
const int pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi::NITERATIONS = 500;

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP1 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity_Terzaghi_TriP1,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity_Terzaghi::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(0, 1), // pressure
            pylith::topology::Field::Discretization(0, 1), // trace strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 1), // porosity
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // fluid_density
            pylith::topology::Field::Discretization(0, 1), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
            pylith::topology::Field::Discretization(0, 1), // biot_coefficient
            pylith::topology::Field::Discretization(0, 1), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 1), // fluid_bulk_modulus

        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 1), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 1), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity_Terzaghi_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP2 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity_Terzaghi_TriP2,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity_Terzaghi::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
            pylith::topology::Field::Discretization(1, 2), // trace strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // fluid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // biot_coefficient
            pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 2), // fluid_bulk_modulus

        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 2), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 2), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity_Terzaghi_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP3 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity_Terzaghi_TriP3,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity_Terzaghi::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
            pylith::topology::Field::Discretization(2, 3), // trace strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 3), // porosity
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // fluid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // biot_coefficient
            pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 3), // fluid_bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 3), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 3), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity_Terzaghi_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ2 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ2,  TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity_Terzaghi::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
            pylith::topology::Field::Discretization(1, 2), // trace strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 2), // porosity
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // fluid_density
            pylith::topology::Field::Discretization(0, 2), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
            pylith::topology::Field::Discretization(0, 2), // biot_coefficient
            pylith::topology::Field::Discretization(0, 2), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 2), // fluid_bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 2), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 2), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ3 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ3,  TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity_Terzaghi::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 3;
        static const pylith::topology::Field::Discretization _solnDiscretizations[3] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
            pylith::topology::Field::Discretization(2, 3), // trace strain
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[9] = {
            pylith::topology::Field::Discretization(0, 3), // porosity
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // fluid_density
            pylith::topology::Field::Discretization(0, 3), // fluid_viscosity
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
            pylith::topology::Field::Discretization(0, 3), // biot_coefficient
            pylith::topology::Field::Discretization(0, 3), // isotropic_permeability
            pylith::topology::Field::Discretization(0, 3), // fluid_bulk_modulus

        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        static const pylith::topology::Field::Discretization _tractionDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 3), // initial_amplitude_tangential
            pylith::topology::Field::Discretization(0, 3), // initial_amplitude_normal
        };
        _data->tractionDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_tractionDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity_Terzaghi_QuadQ3);

// End of file
