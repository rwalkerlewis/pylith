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

#include "PointForceExplosion2D.hh" // Implementation of test data

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include <cmath> // USES M_PI, sin, cos, exp, sqrt

// ------------------------------------------------------------------------------------------------
namespace pylith {
    class _PointForceExplosion2D;
} // pylith

class pylith::_PointForceExplosion2D {
    static pylith::scales::Scales scales;

    // Physical parameters
    static const double WAVELENGTH;
    static const double TIME_SNAPSHOT;

    // Source parameters
    static const double SOURCE_X;
    static const double SOURCE_Y;
    static const double MAGNITUDE;
    static const double FREQUENCY;
    static const double DELAY;
    static const double ORIGIN_TIME;

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

    // Ricker wavelet
    static double rickerWavelet(const double t) {
        const double pi = M_PI;
        const double tau = t - ORIGIN_TIME - DELAY;
        const double pi2_f2_tau2 = (pi * FREQUENCY * tau) * (pi * FREQUENCY * tau);
        return (1.0 - 2.0 * pi2_f2_tau2) * exp(-pi2_f2_tau2);
    } // rickerWavelet

    // Solution subfields - for an explosion source, displacement radiates outward
    // For simplicity, we use a Green's function approximation

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double t) {
        const double pi = M_PI;
        const double velocityScale = scales.getLengthScale() / scales.getTimeScale();
        const double c = vp(x,y) / velocityScale;
        
        // Distance from source
        const double dx = x - SOURCE_X;
        const double dy = y - SOURCE_Y;
        const double r = sqrt(dx*dx + dy*dy);
        
        if (r < 1.0e-10) return 0.0;
        
        // Radial component of displacement from point source
        // u_r ~ M * S(t - r/c) / (4*pi*rho*c^2*r)
        const double rho = density(x, y) / (scales.getRigidityScale() / (scales.getLengthScale() * scales.getLengthScale()));
        const double travelTime = r / c;
        const double retardedTime = t - travelTime;
        
        if (retardedTime < ORIGIN_TIME) return 0.0;
        
        const double stf = rickerWavelet(retardedTime);
        const double amplitude = MAGNITUDE * stf / (4.0 * pi * rho * c * c * r);
        
        return amplitude * dx / r;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double t) {
        const double pi = M_PI;
        const double velocityScale = scales.getLengthScale() / scales.getTimeScale();
        const double c = vp(x,y) / velocityScale;
        
        const double dx = x - SOURCE_X;
        const double dy = y - SOURCE_Y;
        const double r = sqrt(dx*dx + dy*dy);
        
        if (r < 1.0e-10) return 0.0;
        
        const double rho = density(x, y) / (scales.getRigidityScale() / (scales.getLengthScale() * scales.getLengthScale()));
        const double travelTime = r / c;
        const double retardedTime = t - travelTime;
        
        if (retardedTime < ORIGIN_TIME) return 0.0;
        
        const double stf = rickerWavelet(retardedTime);
        const double amplitude = MAGNITUDE * stf / (4.0 * pi * rho * c * c * r);
        
        return amplitude * dy / r;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    // Velocity (time derivative of displacement)
    static double vel_x(const double x,
                        const double y,
                        const double t) {
        // Numerical time derivative
        const double dt = 1.0e-6;
        return (disp_x(x, y, t + dt) - disp_x(x, y, t - dt)) / (2.0 * dt);
    } // vel_x

    static double vel_y(const double x,
                        const double y,
                        const double t) {
        const double dt = 1.0e-6;
        return (disp_y(x, y, t + dt) - disp_y(x, y, t - dt)) / (2.0 * dt);
    } // vel_y

    static const char* vel_units(void) {
        return "m/s";
    } // vel_units

    // Acceleration
    static double acc_x(const double x,
                        const double y,
                        const double t) {
        const double dt = 1.0e-6;
        return (vel_x(x, y, t + dt) - vel_x(x, y, t - dt)) / (2.0 * dt);
    } // acc_x

    static double acc_y(const double x,
                        const double y,
                        const double t) {
        const double dt = 1.0e-6;
        return (vel_y(x, y, t + dt) - vel_y(x, y, t - dt)) / (2.0 * dt);
    } // acc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

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

        s[0] = disp_x(x[0], x[1], t);
        s[1] = disp_y(x[0], x[1], t);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_vel(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = vel_x(x[0], x[1], t);
        s[1] = vel_y(x[0], x[1], t);

        return 0;
    } // solnkernel_vel

    static PetscErrorCode solnkernel_acc(PetscInt spaceDim,
                                         PetscReal t,
                                         const PetscReal x[],
                                         PetscInt numComponents,
                                         PetscScalar* s,
                                         void* context) {
        assert(2 == spaceDim);
        assert(x);
        assert(2 == numComponents);
        assert(s);

        s[0] = acc_x(x[0], x[1], t);
        s[1] = acc_y(x[0], x[1], t);

        return 0;
    } // solnkernel_acc

public:

    static
    TestPointForce_Data* createData(void) {
        TestPointForce_Data* data = new TestPointForce_Data();assert(data);

        data->journalName = "PointForceExplosion2D";
        data->isJacobianLinear = true;
        data->tolerance = 1.0e-4; // Relaxed tolerance due to source regularization

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";
        data->useAsciiMesh = true;

        const double lengthScale = 1.0e+5; // 100 km
        const double velocityScale = vs(0.0, 0.0);
        pylith::scales::ElasticityScales::setDynamicElasticity(&scales, lengthScale, velocityScale);
        data->scales = scales;
        data->formulation = pylith::problems::Physics::DYNAMIC;

        data->t = TIME_SNAPSHOT;
        data->dt = 0.05;

        // Material information
        data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(data->formulation);
        data->material.useBodyForce(false);
        data->rheology.useReferenceState(false);

        data->material.setIdentifier("elasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        // Point source parameters
        const double sourceLocation[2] = {SOURCE_X, SOURCE_Y};
        const double momentTensor[3] = {1.0, 1.0, 0.0}; // Isotropic (explosion)
        data->pointSource.setLocation(sourceLocation, 2);
        data->pointSource.setMomentTensor(momentTensor, 3);
        data->pointSource.setMagnitude(MAGNITUDE);
        data->pointSource.setOriginTime(ORIGIN_TIME);
        data->pointSource.setDominantFrequency(FREQUENCY);
        data->pointSource.setTimeDelay(DELAY);

        // Boundary conditions
        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        pylith::bc::DirichletUserFn* bc = NULL;
        data->bcs.resize(2);
        bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        bc->setUserFnDot(solnkernel_vel);
        data->bcs[0] = bc;

        bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("velocity");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_vel);
        bc->setUserFnDot(solnkernel_acc);
        data->bcs[1] = bc;

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_vel,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        static const pylith::testing::MMSTest::solution_fn _exactSolnDotFns[2] = {
            solnkernel_vel,
            solnkernel_acc,
        };
        data->exactSolnDotFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnDotFns);

        return data;
    } // createData

}; // _PointForceExplosion2D

pylith::scales::Scales pylith::_PointForceExplosion2D::scales;
const double pylith::_PointForceExplosion2D::WAVELENGTH = 1.0e+5;
const double pylith::_PointForceExplosion2D::TIME_SNAPSHOT = 5.0;
const double pylith::_PointForceExplosion2D::SOURCE_X = 0.0;
const double pylith::_PointForceExplosion2D::SOURCE_Y = 0.0;
const double pylith::_PointForceExplosion2D::MAGNITUDE = 1.0e+15;
const double pylith::_PointForceExplosion2D::FREQUENCY = 1.0;
const double pylith::_PointForceExplosion2D::DELAY = 1.0;
const double pylith::_PointForceExplosion2D::ORIGIN_TIME = 0.0;

// ------------------------------------------------------------------------------------------------
pylith::TestPointForce_Data*
pylith::PointForceExplosion2D::TriP1(void) {
    TestPointForce_Data* data = pylith::_PointForceExplosion2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestPointForce_Data*
pylith::PointForceExplosion2D::TriP2(void) {
    TestPointForce_Data* data = pylith::_PointForceExplosion2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestPointForce_Data*
pylith::PointForceExplosion2D::QuadQ1(void) {
    TestPointForce_Data* data = pylith::_PointForceExplosion2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestPointForce_Data*
pylith::PointForceExplosion2D::QuadQ2(void) {
    TestPointForce_Data* data = pylith::_PointForceExplosion2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // vel
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // QuadQ2


// End of file
