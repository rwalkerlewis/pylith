/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/IsotropicLinearPoroelasticityBlackOil.hh"
#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES Isotropic Linear Poroelasticity kernels
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for isotropic, linear poroelasticity plane strain.
// =====================================================================================================================
// ----------------------------------------------------------------------

// ================================= MMS =======================================

// ----------------------------------------------------------------------
// f0u function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_ql_u(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    f0[0] -= (2.0 * shearModulus - biotCoefficient * t);
    f0[1] -= (2.0 * lambda + 4.0 * shearModulus - biotCoefficient * t);
} // f0_quadratic_linear_u


// ----------------------------------------------------------------------
// f0p function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_ql_p(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    PylithScalar sum = 0.0;
    sum += x[0];
    sum += x[1];

    f0[0] -= (sum / biotModulus);
} // f0_quadratic_linear_p


// ----------------------------------------------------------------------
// f0u function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_qt_u(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    f0[0] -= (2.0 * shearModulus - biotCoefficient * PetscCosReal(t));
    f0[1] -= (2.0 * lambda + 4.0 * shearModulus - biotCoefficient * PetscCosReal(t));
} // f0_quadratic_trig_u


// ----------------------------------------------------------------------
// f0p function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_qt_p(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    PylithScalar sum = 0.0;

    sum += x[0];
    sum += x[1];

    f0[0] += PetscSinReal(t)*sum/biotModulus;
} // f0_quadratic_trig_p


// ----------------------------------------------------------------------
// f0u function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_tl_u(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    // f0[0] += PetscSqr(2.0 * PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[0]) * (2.0 * shearModulus + lambda) + 2.0 *
    //          (shearModulus + lambda) - 2.0 * PETSC_PI * biotCoefficient * PetscSinReal(2.0 * PETSC_PI * x[0]) * t;
    // f0[1] += PetscSqr(2.0 * PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[1]) * (2.0 * shearModulus + lambda) - 2.0 *
    //          PETSC_PI * biotCoefficient * PetscSinReal(2.0 * PETSC_PI * x[1]) * t;

    for (PylithInt d = 0; d < _dim-1; ++d) {
        f0[d] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[d])*(2.*shearModulus + lambda) + 2.0*(shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[d])*t;
    }
    f0[_dim-1] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[_dim-1])*(2.*shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[_dim-1])*t;
} // f0_trig_linear_u


// ----------------------------------------------------------------------
// f0p function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0_mms_tl_p(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar kappa = a[aOff[i_isotropicPermeability]] / a[aOff[i_fluidViscosity]];

    PylithScalar sum = 0.0;

    // sum += PetscCosReal(2.0 * PETSC_PI * x[0]) + PetscCosReal(2.0 * PETSC_PI * x[1]);
    for (PylithInt d = 0; d < _dim; ++d) {
        sum += PetscCosReal(2.*PETSC_PI*x[d]);
    }

    f0[0] -= sum / biotModulus - 4.0 * PetscSqr(PETSC_PI) * kappa * sum * t;
} // f0_quadratic_trig_p


// ================================= STD =======================================

// ============================== LHS Residuals ================================

// ----------------------------------------------------------------------
// f0p function for explicit time stepping (dynamic).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_explicit(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += pressure_t/biotModulus;
} // f0p_explicit


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilBlackOilPlaneStrain::f0p_implicit(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += s_t ? (biotCoefficient * trace_strain_t) : 0.0;
    f0[0] += s_t ? (pressure_t / biotModulus) : 0.0;
} // f0p_implicit


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilBlackOilPlaneStrain::f0p_implicit_source(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 4;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_implicit_source_body(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 5;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_body


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_implicit_source_grav(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 5;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_grav


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_implicit_source_grav_body(const PylithInt dim,
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
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 6;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_grav_body


// ----------------------------------------------------------------------
// f1u function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Quasi - Static Case
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1u(const PylithInt dim,
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
                                                                         PylithScalar f1[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar* pressure = &s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff[i_displacement] >= 0);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_displacement] >= 0);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(f1);

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar equivalentPressure = 0;

    for (PylithInt i = 0; i < _phases; i++) {
        equivalentPressure += saturation[i] * pressure[i];
    }

    for (PylithInt c = 0; c < _dim; ++c) {
        for (PylithInt d = 0; d < _dim; ++d) {
            f1[c*_dim+d] -= shearModulus * (displacement_x[c*_dim+d] + displacement_x[d*_dim+c]);
        } // for
        f1[c*_dim+c] -= (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        f1[c*_dim+c] += biotCoefficient*equivalentPressure;
    } // for
} // f1u


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p(const PylithInt dim,
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
                                                                         PylithScalar f1[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOilBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * pressure_x[0];
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * pressure_x[1];

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * pressure_x[0];
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * pressure_x[1];

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * pressure_x[0];
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * pressure_x[1];

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];

} // f1p


// ------------------------------------------*----------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_tensor_permeability(const PylithInt dim,
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
                                                                                             PylithScalar f1[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOilBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim]);
        } // for
    } // for

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];

} // f1p_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_body(const PylithInt dim,
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
                                                                              PylithScalar f1[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_bodyForce = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluidDensity = numA - 13;
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - bodyForce[0]*fluidDensity[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - bodyForce[1]*fluidDensity[0]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[2] - bodyForce[0]*fluidDensity[1]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[3] - bodyForce[1]*fluidDensity[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[4] - bodyForce[0]*fluidDensity[2]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[5] - bodyForce[1]*fluidDensity[2]);

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];

} // f1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH body force, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_body_tensor_permeability(const PylithInt dim,
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
                                                                                                  PylithScalar f1[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_bodyForce = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluidDensity = numA - 13;
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - bodyForce[j]*fluidDensity[0]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - bodyForce[j]*fluidDensity[1]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - bodyForce[j]*fluidDensity[2]);
        } // for
    } // for

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];
} // f1p_gravity_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_gravity(const PylithInt dim,
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
                                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
    } // for

} // f1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                                     PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *gravityField = &a[aOff[i_gravityField]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
        } // for
    } // for

} // f1p_gravity_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_body_gravity(const PylithInt dim,
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
                                                                                      PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluid_density = 1;
    const PylithInt i_fluid_viscosity = 2;

    const PylithInt i_body_force = 4;
    const PylithInt i_gravity_field = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropic_permeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_fluid_viscosity] >= 0);
    assert(aOff[i_isotropic_permeability] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluid_density = a[aOff[i_fluid_density]];
    const PylithScalar fluid_viscosity = a[aOff[i_fluid_viscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

    const PylithScalar isotropic_permeability = a[aOff[i_isotropic_permeability]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropic_permeability / fluid_viscosity) * (pressure_x[d] - body_force[d] - fluid_density * gravity_field[d]);
    } // for

} // f1p_body_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1p_body_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                                          PylithScalar f1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;

    const PylithInt i_body_force = 4;
    const PylithInt i_gravity_field = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j] - fluidDensity * gravity_field[j]);
        } // for
    } // for

} // f1p_body_gravity_tensor_permeability


// =============================== LHS Jacobian ================================

// ---------------------------------------------------------------------------------------------------------------------
/* Jf3_uu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf3uu(const PylithInt dim,
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
                                                                           const PylithReal s_tshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf3[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.

    // Incoming auxiliary fields.

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(Jf3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            Jf3[((i*_dim + i)*_dim + j)*_dim + j] -= shearModulus;
            Jf3[((i*_dim + j)*_dim + j)*_dim + i] -= shearModulus;
        }
    }

} // Jf3uu


// -----------------------------------------------------------------------------
// Jf2up function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf2up(const PylithInt dim,
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
                                                                           const PylithReal utshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf2[]) {
    const PylithInt _dim = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf2);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf2[d*_dim+d] += biotCoefficient;
    } // for
} // Jf2up


// -----------------------------------------------------------------------------
// Jf2ue function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf2ue(const PylithInt dim,
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
                                                                           const PylithReal utshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf2[]) {
    const PylithInt _dim = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_drainedBulkModulus] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(Jf2);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf2[d*_dim+d] -= lambda;
    } // for
} // Jf2ue


// ----------------------------------------------------------------------
/* Jf3pp entry function for isotropic linear poroelasticity. */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf3pp(const PylithInt dim,
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
                                                                           const PylithReal utshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf3[]) {
    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(Jf3);

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf3[d*_dim+d] += isotropicPermeablity / fluidViscosity;
    }
} // Jf3pp


// ----------------------------------------------------------------------
/* Jf3pp entry function for isotropic linear poroelasticity, permeability in tensor form */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf3pp_tensor_permeability(const PylithInt dim,
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
                                                                                               const PylithReal utshift,
                                                                                               const PylithScalar x[],
                                                                                               const PylithInt numConstants,
                                                                                               const PylithScalar constants[],
                                                                                               PylithScalar Jf3[]) {
    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(Jf3);

    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            Jf3[i*_dim+j] += tensorPermeability[i*_dim+j]/fluidViscosity;
        } // for
    } // for
} // Jf3pp_tensorPermeability


// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0pp(const PylithInt dim,
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
                                                                           const PylithReal utshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(Jf0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] += utshift / biotModulus;
} // Jf0pp


// -----------------------------------------------------------------------------
// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0pe(const PylithInt dim,
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
                                                                           const PylithReal utshift,
                                                                           const PylithScalar x[],
                                                                           const PylithInt numConstants,
                                                                           const PylithScalar constants[],
                                                                           PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticity
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] += utshift * biotCoefficient;
} // Jf0pe


// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0ppdot(const PylithInt dim,
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
                                                                              const PylithReal utshift,
                                                                              const PylithScalar x[],
                                                                              const PylithInt numConstants,
                                                                              const PylithScalar constants[],
                                                                              PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(Jf0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] += 1.0 / biotModulus;
} // Jf0ppdot


// -----------------------------------------------------------------------------
// Jf0pedot function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0pedot(const PylithInt dim,
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
                                                                              const PylithReal utshift,
                                                                              const PylithScalar x[],
                                                                              const PylithInt numConstants,
                                                                              const PylithScalar constants[],
                                                                              PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] += biotCoefficient;
} // Jf0pedot


// ============================== RHS Residual =================================

// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g0p(const PylithInt dim,
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
                                                                         PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_implicit


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g0p_source(const PylithInt dim,
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
                                                                                PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 4;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g0p_source_body(const PylithInt dim,
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
                                                                                     PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_body


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g0p_source_grav(const PylithInt dim,
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
                                                                                     PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_grav


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g0p_source_grav_body(const PylithInt dim,
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
                                                                                          PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 6;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_grav_body


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1p_gravity(const PylithInt dim,
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
                                                                                 PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1[d] -= (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity*gravityField[d]);
    } // for

} // g1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1p_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                                     PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1[i] -= (tensorPermeability[i*_dim+j] / fluidViscosity) * (pressure_x[j] - fluidDensity*gravityField[j]);
        } // for
    } // for

} // g1p_gravity_tensor_permeability


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1p(const PylithInt dim,
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
                                                                         PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1[d] -= (isotropicPermeability / fluidViscosity) * pressure_x[d];
    } // for
} // g1p


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_tensor_permeability(const PylithInt dim,
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
                                                                                     PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPermeability[4];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1[i] -= (tensorPermeability[i*_dim+j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // g1p_tensor_permeability


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Dynamic Case
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1v(const PylithInt dim,
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
                                                                         PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain += displacement_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt c = 0; c < _dim; ++c) {
        for (PylithInt d = 0; d < _dim; ++d) {
            g1[c*dim+d] -= shearModulus * (displacement_x[c*_dim+d] + displacement_x[d*_dim+c]);
        } // for
        g1[c*dim+c] -= (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        g1[c*dim+c] += biotCoefficient*pressure;
    } // for
} // g1v


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1v_refstate(const PylithInt dim,
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
                                                                                  PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain += displacement_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
            refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
        } // for
    } // for

    for (PylithInt i = 0; i < _dim; ++i) {
        g1[i*_dim+i] -= (meanStress - alphaPres);
        g1[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            g1[i*_dim+j] -= shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for
} // g1v_refstate


// RHS Jacobians

// ========================== Helper Kernels ===================================

// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [displacement(dim), ...]
 * Auxiliary fields: [density(1), ..., shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::cauchyStress(const PylithInt dim,
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
                                                                                  PylithScalar stressVector[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    // Create and populate stress tensor

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] += shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]);
        } // for
        stressTensor[i*_dim+i] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        stressTensor[i*_dim+i] -= biotCoefficient*pressure;
    } // for

    // Construct stress vector
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0] + stressTensor[1*_dim+1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [displacement(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::cauchyStress_refstate(const PylithInt dim,
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
                                                                                           PylithScalar stressVector[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
            refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
        } // for
    } // for

    // Create and populate stress tensor

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] -= (meanStress - alphaPres);
        stressTensor[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] -= shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for

    // Generate stress vector

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = refStress[2] +  0.5*lambda/(lambda+shearModulus) *
                                   (stressTensor[0*_dim+0]-refStress[0] + stressTensor[1*_dim+1]-refStress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate


// ========================== Update Kernels ===================================

// ---------------------------------------------------------------------------------------------------------------------
/* Update porosity for a linear poroelastic material, implicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::updatePorosityImplicit(const PylithInt dim,
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
                                                                                            PylithScalar porosity[]) {
    const PylithInt _dim = 2;
    // Incoming solution fields.
    const PylithInt i_pressure_t = 4;
    const PylithInt i_trace_strain_t = 5;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosityPrev = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosityPrev] >= 0);
    assert(porosity);

    // Do stuff
    const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
    const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

    // Update porosity
    porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                       ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                       drainedBulkModulus * pressure_t);

} // updatePorosityImplicit


// ---------------------------------------------------------------------------------------------------------------------
/* Update porosity for a linear poroelastic material, explicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::updatePorosityExplicit(const PylithInt dim,
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
                                                                                            PylithScalar porosity[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosity] >= 0);
    assert(porosity);

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

    // Do stuff
    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * pressure_t);
} // updatePorosityExplicit


// =====================================================================================================================
// Kernels for isotropic, linear poroelasticity in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------

// ================================= MMS =======================================

// ----------------------------------------------------------------------
// f0u function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_ql_u(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_drainedBulkModulus] >= 0);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(f0);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    for (PylithInt d = 0; d < _dim-1; ++d) {
        f0[d] -= 2.0*shearModulus - biotCoefficient*t;
    }
    f0[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*t;
} // f0_quadratic_linear_u


// ----------------------------------------------------------------------
// f0p function for quadratic space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_ql_p(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.

    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    PylithScalar sum = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        sum += x[d];
    }
    f0[0] -= sum/biotModulus;
} // f0_quadratic_linear_p


// ----------------------------------------------------------------------
// f0u function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_qt_u(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    for (PylithInt d = 0; d < _dim-1; ++d) {
        f0[d] -= 2.0*shearModulus - biotCoefficient*PetscCosReal(t);
    }
    f0[_dim-1] -= 2.0*lambda + 4.0*shearModulus - biotCoefficient*PetscCosReal(t);
} // f0_quadratic_trig_u


// ----------------------------------------------------------------------
// f0p function for quadratic space and trigonometric time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_qt_p(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.

    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    PylithScalar sum = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        sum += x[d];
    }

    f0[0] += PetscSinReal(t)*sum/biotModulus;
} // f0_quadratic_trig_p


// ----------------------------------------------------------------------
// f0u function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_tl_u(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;

    for (PylithInt d = 0; d < _dim-1; ++d) {
        f0[d] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[d])*(2.*shearModulus + lambda) + 2.0*(shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[d])*t;
    }
    f0[_dim-1] += PetscSqr(2.*PETSC_PI)*PetscSinReal(2.*PETSC_PI*x[_dim-1])*(2.*shearModulus + lambda) - 2.*PETSC_PI*biotCoefficient*PetscSinReal(2.*PETSC_PI*x[_dim-1])*t;
} // f0_trig_linear_u


// ----------------------------------------------------------------------
// f0p function for trigonometric space and linear time MMS.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0_mms_tl_p(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar kappa = a[aOff[i_isotropicPermeability]] / a[aOff[i_fluidViscosity]];
    PylithScalar sum = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        sum += PetscCosReal(2.*PETSC_PI*x[d]);
    }

    f0[0] -= sum/biotModulus - 4*PetscSqr(PETSC_PI)*kappa*sum*t;
} // f0_quadratic_trig_p


// ================================= STD =======================================

// ================================= LHS =======================================

// ----------------------------------------------------------------------
// f0p function for explicit time stepping (dynamic).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_explicit(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;

    // Incoming re-packed auxiliary field.

    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += pressure_t/biotModulus;
} // f0p_explicit


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_implicit(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
} // f0p_implicit


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_implicit_source(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 4;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_implicit_source_body(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 5;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_body


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_implicit_source_grav(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 5;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_grav


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0p_implicit_source_grav_body(const PylithInt dim,
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
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;
    const PylithInt i_source = 6;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(sOff_x[i_trace_strain] >= 0);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(aOff[i_biotModulus] >= 0);
    assert(aOff[i_source] >= 0);
    assert(f0);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    f0[0] += biotCoefficient*trace_strain_t;
    f0[0] += pressure_t/biotModulus;
    f0[0] -= source;
} // f0p_implicit_source_grav_body


// ----------------------------------------------------------------------
// f1u function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Quasi - Static Case
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1u(const PylithInt dim,
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
                                                                PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_displacement] >= 0);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(f1);

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt c = 0; c < _dim; ++c) {
        for (PylithInt d = 0; d < _dim; ++d) {
            f1[c*_dim+d] -= shearModulus * (displacement_x[c*_dim+d] + displacement_x[d*_dim+c]);
        } // for
        f1[c*_dim+c] -= (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        f1[c*_dim+c] += biotCoefficient*pressure;
    } // for
} // f1u


// ----------------------------------------------------------------------
// f1u function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1u_refstate(const PylithInt dim,
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
                                                                         PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanrstress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanrstress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
            refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
        } // for
    } // for

    for (PylithInt i = 0; i < _dim; ++i) {
        f1[i*_dim+i] -= (meanStress - alphaPres);
        f1[i*_dim+i] -= refStress[i*_dim+i] - meanrstress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            f1[i*_dim+j] -= shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]) - refStrain[i*_dim+j];
        } // for
    } // for
} // f1u_refstate


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p(const PylithInt dim,
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
                                                                PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d]);
    } // for
} // f1p


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_tensor_permeability(const PylithInt dim,
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
                                                                                    PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // f1p_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_body(const PylithInt dim,
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
                                                                     PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;

    const PylithInt i_body_force = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - body_force[d]);
    } // for

} // f1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH body force, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_body_tensor_permeability(const PylithInt dim,
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
                                                                                         PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_body_force = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j]);
        } // for
    } // for

} // f1p_gravity_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_gravity(const PylithInt dim,
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
                                                                        PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
    } // for

} // f1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                            PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *gravityField = &a[aOff[i_gravityField]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
        } // for
    } // for

} // f1p_gravity_tensor_permeability


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_body_gravity(const PylithInt dim,
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
                                                                             PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluid_density = 1;
    const PylithInt i_fluid_viscosity = 2;

    const PylithInt i_body_force = 4;
    const PylithInt i_gravity_field = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropic_permeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_fluid_viscosity] >= 0);
    assert(aOff[i_isotropic_permeability] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluid_density = a[aOff[i_fluid_density]];
    const PylithScalar fluid_viscosity = a[aOff[i_fluid_viscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

    const PylithScalar isotropic_permeability = a[aOff[i_isotropic_permeability]];

    for (PylithInt d = 0; d < _dim; ++d) {
        f1[d] += (isotropic_permeability / fluid_viscosity) * (pressure_x[d] - body_force[d] - fluid_density * gravity_field[d]);
    } // for

} // f1p_body_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f1p_body_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                                 PylithScalar f1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;

    const PylithInt i_body_force = 4;
    const PylithInt i_gravity_field = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 4);
    assert(sOff_x);
    assert(sOff_x[i_pressure] >= 0);
    assert(aOff);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    const PylithScalar *body_force = &a[aOff[i_body_force]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

    const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j] - fluidDensity * gravity_field[j]);
        } // for
    } // for

} // f1p_body_gravity_tensor_permeability


// =============================== LHS Jacobian ================================

// ---------------------------------------------------------------------------------------------------------------------
/* Jg3_vu entry function for isotropic linear poroelasticity WITHOUT reference stress and reference strain.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22, C1111=lambda+2mu, C1122=lambda.
 *
 * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl
 *   + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf3uu(const PylithInt dim,
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
                                                                  const PylithReal s_tshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf3[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.

    // Incoming auxiliary fields.

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 1);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_shearModulus] >= 0);
    assert(Jf3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            Jf3[((i*_dim + i)*_dim + j)*_dim + j] -= shearModulus;
            Jf3[((i*_dim + j)*_dim + j)*_dim + i] -= shearModulus;
        }
    }

} // Jf3uu


// -----------------------------------------------------------------------------
// Jf2up function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf2up(const PylithInt dim,
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
                                                                  const PylithReal utshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf2[]) {
    const PylithInt _dim = 3;

    // Isotropic Linear Poroelasticity
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf2);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf2[d*_dim+d] += biotCoefficient;
    } // for
} // Jf2up


// -----------------------------------------------------------------------------
// Jf2ue function for isotropic linear poroelasticity.

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf2ue(const PylithInt dim,
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
                                                                  const PylithReal utshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf2[]) {
    const PylithInt _dim = 3;

    // Isotropic Linear Poroelasticity
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_drainedBulkModulus] >= 0);
    assert(aOff[i_shearModulus] >= 0);
    assert(Jf2);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf2[d*_dim+d] -= drainedBulkModulus - (2.0*shearModulus) / 3.0;
    } // for
} // Jf2ue


// ----------------------------------------------------------------------
/* Jf3pp entry function for isotropic linear poroelasticity. */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf3pp(const PylithInt dim,
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
                                                                  const PylithReal utshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf3[]) {
    const PylithInt _dim = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(Jf3);

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf3[d*_dim+d] += isotropicPermeablity/fluidViscosity;
    }
} // Jf3pp


// ----------------------------------------------------------------------
/* Jf3pp entry function for isotropic linear poroelasticity, permeability in tensor form */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf3pp_tensor_permeability(const PylithInt dim,
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
                                                                                      const PylithReal utshift,
                                                                                      const PylithScalar x[],
                                                                                      const PylithInt numConstants,
                                                                                      const PylithScalar constants[],
                                                                                      PylithScalar Jf3[]) {
    const PylithInt _dim = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // Isotropic Linear Poroelasticity
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(Jf3);

    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            Jf3[i*_dim+j] += tensorPermeability[i*_dim+j]/fluidViscosity;
        } // for
    } // for
} // Jf3pp_tensorPermeability


// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0pp(const PylithInt dim,
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
                                                                  const PylithReal utshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf0[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(Jf0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] = utshift / biotModulus;
} // Jf0pp


// -----------------------------------------------------------------------------
// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0pe(const PylithInt dim,
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
                                                                  const PylithReal utshift,
                                                                  const PylithScalar x[],
                                                                  const PylithInt numConstants,
                                                                  const PylithScalar constants[],
                                                                  PylithScalar Jf0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] += utshift * biotCoefficient;
} // Jf0pe


// -----------------------------------------------------------------------------
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0ppdot(const PylithInt dim,
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
                                                                     const PylithReal utshift,
                                                                     const PylithScalar x[],
                                                                     const PylithInt numConstants,
                                                                     const PylithScalar constants[],
                                                                     PylithScalar Jf0[]) {
    const PylithInt _dim = 3;

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotModulus] >= 0);
    assert(Jf0);

    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    Jf0[0] += 1.0 / biotModulus;
} // Jf0ppdot


// -----------------------------------------------------------------------------
// Jf0pedot function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0pedot(const PylithInt dim,
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
                                                                     const PylithReal utshift,
                                                                     const PylithScalar x[],
                                                                     const PylithInt numConstants,
                                                                     const PylithScalar constants[],
                                                                     PylithScalar Jf0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    Jf0[0] += biotCoefficient;
} // Jf0pedot


// ============================== RHS Residual =================================

// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g0p(const PylithInt dim,
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
                                                                PylithScalar g0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_implicit


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g0p_source(const PylithInt dim,
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
                                                                       PylithScalar g0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 4;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g0p_source_body(const PylithInt dim,
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
                                                                            PylithScalar g0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_body


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g0p_source_grav(const PylithInt dim,
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
                                                                            PylithScalar g0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_grav


// ----------------------------------------------------------------------
// g0p function for generic poroelasticity terms.
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g0p_source_grav_body(const PylithInt dim,
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
                                                                                 PylithScalar g0[]) {
    const PylithInt _dim = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 6;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    g0[0] += source;
    g0[0] -= biotCoefficient*trace_strain_t;
} // g0p_source_grav_body


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1p_gravity(const PylithInt dim,
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
                                                                        PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1[d] -= (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity*gravityField[d]);
    } // for

} // g1p_gravity


// -----------------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITH gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1p_gravity_tensor_permeability(const PylithInt dim,
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
                                                                                            PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1[i] -= (tensorPermeability[i*_dim+j] / fluidViscosity) * (pressure_x[j] - fluidDensity*gravityField[j]);
        } // for
    } // for

} // g1p_gravity_tensor_permeability


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1p(const PylithInt dim,
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
                                                                PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_isotropicPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    for (PylithInt d = 0; d < _dim; ++d) {
        g1[d] -= (isotropicPermeability / fluidViscosity) * pressure_x[d];
    } // for
} // g1p


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for isotropic linear
 * poroelasticity WITHOUT gravity, permeability represented in tensor form.
 *
 * darcyFlow = -k / mu_f (det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1p_tensor_permeability(const PylithInt dim,
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
                                                                                    PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_tensorPermeability = numA - 1;

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

    PylithScalar tensorPermeability[9];
    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            g1[i] -= (tensorPermeability[i*_dim+j] / fluidViscosity) * (pressure_x[j]);
        } // for
    } // for
} // g1p_tensor_permeability


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity WITHOUT reference stress and strain.
// Dynamic Case
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1v(const PylithInt dim,
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
                                                                PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain += displacement_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    for (PylithInt c = 0; c < _dim; ++c) {
        for (PylithInt d = 0; d < _dim; ++d) {
            g1[c*dim+d] -= shearModulus * (displacement_x[c*_dim+d] + displacement_x[d*_dim+c]);
        } // for
        g1[c*dim+c] -= (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        g1[c*dim+c] += biotCoefficient*pressure;
    } // for
} // g1v


// ----------------------------------------------------------------------
// g1v function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::g1v_refstate(const PylithInt dim,
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
                                                                         PylithScalar g1[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];

    PylithScalar trace_strain = 0.0;
    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain += displacement_x[d*_dim+d];
    }

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
            refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
        } // for
    } // for

    for (PylithInt i = 0; i < _dim; ++i) {
        g1[i*_dim+i] -= (meanStress - alphaPres);
        g1[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            g1[i*_dim+j] -= shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for
} // g1v_refstate


// ========================== Helper Kernels ===================================

// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [displacement(dim), ...]
 * Auxiliary fields: [density(1), ..., shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::cauchyStress(const PylithInt dim,
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
                                                                         PylithScalar stressVector[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];

    // Create and populate stress tensor

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] += shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]);
        } // for
        stressTensor[i*_dim+i] += (drainedBulkModulus - (2.0*shearModulus)/3.0) * trace_strain;
        stressTensor[i*_dim+i] -= biotCoefficient*pressure;
    } // for

    // Construct stress vector
    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0*_dim+0] + stressTensor[1*_dim+1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress


// ---------------------------------------------------------------------------------------------------------------------
/* Calculate stress vector for isotropic linear
 * poroelasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [displacement(dim), ...]
 * Auxiliary fields: [density(1), ..., refstress(4), refstrain(4), shear_modulus(1), bulk_modulus(1)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::cauchyStress_refstate(const PylithInt dim,
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
                                                                                  PylithScalar stressVector[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_rstress = numA - 7;
    const PylithInt i_rstrain = numA - 6;
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_biotModulus = numA - 2;

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar biotModulus = a[aOff[i_biotModulus]];
    const PylithScalar* refStress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refStrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithScalar ref_trace_strain = refStrain[0] + refStrain[1] + refStrain[2];
    const PylithScalar meanRefStress = (refStress[0] + refStress[1] + refStress[2]) / 3.0;
    const PylithScalar meanStress = meanRefStress + drainedBulkModulus * (trace_strain - ref_trace_strain);
    const PylithScalar alphaPres = biotCoefficient * pressure;
    const PylithReal traceTerm = (-2.0/3.0)*shearModulus * trace_strain;

    // Convert reference vectors to refrence tensors
    PylithScalar refStressTensor[_dim*_dim];
    PylithScalar refStrainTensor[_dim*_dim];
    PylithInt refTensorPos[9] = {0, 3, 5, 3, 1, 4, 5, 4, 2};

    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; ++j) {
            refStressTensor[i*_dim+j] = refStress[refTensorPos[i*_dim+j]];
            refStrainTensor[i*_dim+j] = refStrain[refTensorPos[i*_dim+j]];
        } // for
    } // for

    // Create and populate stress tensor

    PylithScalar stressTensor[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for (PylithInt i = 0; i < _dim; ++i) {
        stressTensor[i*_dim+i] -= (meanStress - alphaPres);
        stressTensor[i*_dim+i] -= refStressTensor[i*_dim+i] - meanRefStress + traceTerm;
        for (PylithInt j = 0; j < _dim; ++j) {
            stressTensor[i*_dim+j] -= shearModulus * (displacement_x[i*_dim+j] + displacement_x[j*_dim+i]) - refStrainTensor[i*_dim+j];
        } // for
    } // for

    // Generate stress vector

    const PylithScalar lambda = drainedBulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = refStress[2] +  0.5*lambda/(lambda+shearModulus) *
                                   (stressTensor[0*_dim+0]-refStress[0] + stressTensor[1*_dim+1]-refStress[1]);

    stressVector[0] = stressTensor[0*_dim+0]; // stress_xx
    stressVector[1] = stressTensor[1*_dim+1]; // stress_yy
    stressVector[2] = stress_zz;
    stressVector[3] = stressTensor[0*_dim+1]; // stress_xy
} // cauchyStress_refstate


// ========================== Update Kernels ===================================

// ---------------------------------------------------------------------------------------------------------------------
/* Update porosity for a linear poroelastic material, implicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::updatePorosityImplicit(const PylithInt dim,
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
                                                                                   PylithScalar porosity[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_pressure_t = 4;
    const PylithInt i_trace_strain_t = 5;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosity] >= 0);
    assert(porosity);

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

    // Do stuff
    const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
    const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * pressure_t);
} // updatePorosityImplicit


// ---------------------------------------------------------------------------------------------------------------------
/* Update porosity for a linear poroelastic material, explicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::updatePorosityExplicit(const PylithInt dim,
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
                                                                                   PylithScalar porosity[]) {
    const PylithInt _dim = 3;

    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosity] >= 0);
    assert(porosity);

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

    // Do stuff
    const PylithScalar pressure_t = s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * pressure_t);
} // updatePorosityExplicit


// End of file
