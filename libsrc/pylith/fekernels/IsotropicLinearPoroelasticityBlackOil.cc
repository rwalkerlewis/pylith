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

// ================================= STD =======================================

// ============================== LHS Residuals ================================

// ---------------------------------------------------------------------------------------------------------------------
// f0u placeholder function for multiphaseporoelasticity equation
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0u_implicit(const PylithInt dim,
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
    // Incoming solution fields.

    // Incoming auxiliary fields.

    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += 0.0;
        // PetscPrintf(PETSC_COMM_WORLD, "f0u[%i]: %f\n",i, f0[i]);
    } // for
} // f0u

// ---------------------------------------------------------------------------------------------------------------------
// f0u placeholder function for multiphaseporoelasticity equation
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0u_body_implicit(const PylithInt dim,
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
    // Incoming solution fields.

    // Incoming auxiliary fields

    // MultiphasePoroelasticity

    // 3 + n
    const PylithInt i_bodyForce = 4;

    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += bodyForce[i];
    } // for
} // f0u_body_implicit

// ---------------------------------------------------------------------------------------------------------------------
// f0u_grav_implicit - f0 function for generic multiphaseporoelasticity terms ( + grav body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0u_grav_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;                                                          
    // Incoming auxililary fields.

    // MultiphasePoroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_porosity = 1;

    // 3 + n
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluid_density = numA - 13;
    const PylithInt i_threePhaseSaturation = numA - 7;

    const PylithScalar* fluidDensity = &a[aOff[i_fluid_density]];
    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar solidDensity = a[aOff[i_solid_density]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    PylithScalar equivalentFluidDensity = 0.0, bulkDensity = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        equivalentFluidDensity += saturation[i]*fluidDensity[i];
    }

    bulkDensity = (1.0 - porosity)*solidDensity + porosity*equivalentFluidDensity;

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += bulkDensity * gravityField[i];
    } // for

} // f0u_grav_implicit

// ---------------------------------------------------------------------------------------------------------------------
//f0u_body_grav_implicit - f0 function for generic multiphaseporoelasticity terms ( + grav body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0u_body_grav_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;                                                          
    // Incoming auxililary fields.

    // MultiphasePoroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_porosity = 1;

    // 3 + n
    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluid_density = numA - 13;
    const PylithInt i_threePhaseSaturation = numA - 7;

    const PylithScalar* fluidDensity = &a[aOff[i_fluid_density]];
    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar solidDensity = a[aOff[i_solid_density]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    PylithScalar equivalentFluidDensity = 0.0, bulkDensity = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        equivalentFluidDensity += saturation[i]*fluidDensity[i];
    }

    bulkDensity = (1.0 - porosity)*solidDensity + porosity*equivalentFluidDensity;

    // gravity field
    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += bulkDensity * gravityField[i];
    } // for

    // body force
    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += bodyForce[i];
    } // for    

} // f0u_body_grav_implicit


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
    //const PylithInt _dim = 2;

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
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t;
    } // for
} // f0p_implicit


// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (body forces).
// \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
// \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}

void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0p_implicit_source(const PylithInt dim,
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for

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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 6;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
} // f0p_implicit_source_grav_body

// ---------------------------------------------------------------------------------------------------------------------
// f0pdot function for multiphaseporoelasticity equation, implicit time stepping, quasistatic.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f0pdot_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;
    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_pdot = 4;

    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_pdot] >= 0);
    assert(s);
    assert(s_t);

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]]; // pressure_t
    const PylithScalar* pdot = &s[sOff[i_pdot]]; // pdot

    for (PylithInt i = 0; i < _phases; ++i) {
        f0[i] += pressure_t[i];
        f0[i] -= pdot[i];
    }
} // f0pdot_implicit

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
// f1u function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::f1u_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar pressure = s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    // Incoming auxiliary fields.

    // IsotropicLinearPoroelasticity
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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
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
    //const PylithInt _phases = 3;

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_bodyForce = 4;

    // IsotropicLinearPoroelasticityBlackOil
    // const PylithInt i_fluidDensity = numA - 13;
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

    //const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - bodyForce[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - bodyForce[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[2] - bodyForce[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[3] - bodyForce[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[4] - bodyForce[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[5] - bodyForce[1]);

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
    //const PylithInt _phases = 3;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    //const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
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
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - bodyForce[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - bodyForce[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - bodyForce[j]);
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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOilBlackOil
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
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - fluidDensity[0] * gravityField[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - fluidDensity[0] * gravityField[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[0] - fluidDensity[1] * gravityField[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[1] - fluidDensity[1] * gravityField[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[0] - fluidDensity[2] * gravityField[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[1] - fluidDensity[2] * gravityField[1]);

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_gravityField = 4;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];
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
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - fluidDensity[0] * gravityField[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - fluidDensity[1] * gravityField[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - fluidDensity[2] * gravityField[j]);
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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_threePhaseViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - bodyForce[0] - fluidDensity[0] * gravityField[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - bodyForce[1] - fluidDensity[0] * gravityField[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[0] - bodyForce[0] - fluidDensity[1] * gravityField[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[1] - bodyForce[1] - fluidDensity[1] * gravityField[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[0] - bodyForce[0] - fluidDensity[2] * gravityField[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[1] - bodyForce[1] - fluidDensity[2] * gravityField[1]);

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
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
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - bodyForce[j] - fluidDensity[0] * gravityField[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - bodyForce[j] - fluidDensity[1] * gravityField[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - bodyForce[j] - fluidDensity[2] * gravityField[j]);
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
    const PylithInt _phases = 3;

    // Isotropic Linear Poroelasticity

    const PylithInt i_saturation = numA - 7;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf2);

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar equivalentBiotCoefficient = 0;

    for (PylithInt i = 0; i < _phases; i++) {
        equivalentBiotCoefficient += saturation[i] * biotCoefficient;
    }

    for (PylithInt d = 0; d < _dim; ++d) {
        Jf2[d*_dim+d] += equivalentBiotCoefficient;
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
    const PylithInt _phases = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(Jf3);

    const PylithScalar *formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar *relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar *fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];

    PylithScalar solutionRatios[_dim*_dim], coeff[_dim*_dim];
    PylithScalar tensorPermeability[_dim*_dim];

    tensorPermeability[0] = 0.0;
    tensorPermeability[1] = 0.0;
    tensorPermeability[2] = 0.0;
    tensorPermeability[3] = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        tensorPermeability[d*_dim+d] += isotropicPermeablity;
    }

    solutionRatios[0] = 1.0; // R_ww
    solutionRatios[1] = 0.0; // R_wo
    solutionRatios[2] = 0.0; // R_wg
    solutionRatios[3] = 0.0; // R_ow
    solutionRatios[4] = 1.0; // R_oo
    solutionRatios[5] = a[aOff[i_solutionOilGasRatio]]; // R_og
    solutionRatios[6] = 0.0; // R_gw
    solutionRatios[7] = a[aOff[i_solutionGasOilRatio]]; // R_go
    solutionRatios[8] = 1.0; // R_gg

    coeff[0] = (solutionRatios[0] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ww
    coeff[1] = (solutionRatios[1] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // wo
    coeff[2] = (solutionRatios[2] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // wg

    coeff[3] = (solutionRatios[3] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ow
    coeff[4] = (solutionRatios[4] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // oo
    coeff[5] = (solutionRatios[5] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // og

    coeff[6] = (solutionRatios[6] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // gw
    coeff[7] = (solutionRatios[7] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // go
    coeff[8] = (solutionRatios[8] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // gg

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            for (PylithInt k = 0; k < _dim; ++k) {
                for (PylithInt l = 0; l < _dim; ++l) {
                    Jf3[((i * _phases + j) * _dim + k) * _dim + l] -= coeff[i*_phases+j] * tensorPermeability[k*_dim+l];
                }
            }
        }
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
    const PylithInt _phases = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(Jf3);

    const PylithScalar *formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar *relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar *fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar solutionRatios[_dim*_dim], coeff[_dim*_dim];
    PylithScalar tensorPermeability[_dim*_dim];

    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[1];
    tensorPermeability[3] = vectorPermeability[3];

    solutionRatios[0] = 1.0; // R_ww
    solutionRatios[1] = 0.0; // R_wo
    solutionRatios[2] = 0.0; // R_wg
    solutionRatios[3] = 0.0; // R_ow
    solutionRatios[4] = 1.0; // R_oo
    solutionRatios[5] = a[aOff[i_solutionOilGasRatio]]; // R_og
    solutionRatios[6] = 0.0; // R_gw
    solutionRatios[7] = a[aOff[i_solutionGasOilRatio]]; // R_go
    solutionRatios[8] = 1.0; // R_gg

    coeff[0] = (solutionRatios[0] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ww
    coeff[1] = (solutionRatios[1] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // wo
    coeff[2] = (solutionRatios[2] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // wg

    coeff[3] = (solutionRatios[3] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ow
    coeff[4] = (solutionRatios[4] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // oo
    coeff[5] = (solutionRatios[5] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // og

    coeff[6] = (solutionRatios[6] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // gw
    coeff[7] = (solutionRatios[7] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // go
    coeff[8] = (solutionRatios[8] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // gg

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            for (PylithInt k = 0; k < _dim; ++k) {
                for (PylithInt l = 0; l < _dim; ++l) {
                    Jf3[((i * _phases + j) * _dim + k) * _dim + l] -= coeff[i*_phases+j] * tensorPermeability[k*_dim+l];
                }
            }
        }
    }
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
    const PylithInt _phases = 3;

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Jf0[i*_phases + j] = utshift * N_tensor[i*_phases + j];
        }
    }
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
    const PylithInt _phases = 3;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    // Jf0pe(n_p, n_e)
    for (PylithInt i = 0; i < _phases; ++i) {
        Jf0[i] += utshift * saturation[i] * biotCoefficient;
    }
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
    const PylithInt _phases = 3;

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Jf0[i*_phases + j] = N_tensor[i*_phases + j];
        }
    }
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
    const PylithInt _phases = 3;

    // Incoming re-packed auxiliary field.
    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_biotCoefficient = numA - 3;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    // Jf0pe(n_p, n_e)
    for (PylithInt i = 0; i < _phases; ++i) {
        Jf0[i] += saturation[i] * biotCoefficient;
    }
} // Jf0pedot


// -----------------------------------------------------------------------------
// Jf0pdotp - Jf0 function for isotropic linear multiphaseporoelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0pdotp(const PylithInt dim,
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
                                                                              PylithScalar Jf0[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    assert(aOff);
    assert(a);

    for (PylithInt p = 0; p < _phases; ++p) {
        Jf0[p*_dim+p] -= s_tshift;
    } // for
} // Jg0pdotp


// -----------------------------------------------------------------------------
// Jf0pdotpdot - Jf0 function for isotropic linear multiphaseporoelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::Jf0pdotpdot(const PylithInt dim,
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
                                                                                 PylithScalar Jf0[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;

    assert(aOff);
    assert(a);

    for (PylithInt p = 0; p < _phases; ++p) {
        Jf0[p*_dim+p] -= s_tshift;
    } // for
} // Jg0pdotpdot


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
    //const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;

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
    //const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 4;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 6;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::g1p_tensor_permeability(const PylithInt dim,
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


    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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
/* Update porosity for a linear poroelastic material, implicit. ASSUMES NO CAPILLARY PRESSURE
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
    const PylithInt _phases = 3;
    // Incoming solution fields.
    const PylithInt i_pressure_t = 4;
    const PylithInt i_trace_strain_t = 5;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosityPrev = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_threePhaseSaturation = numA - 7;
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
    const PylithScalar* pressure_t = &s[sOff[i_pressure_t]];
    const PylithScalar trace_strain_t = s[sOff[i_trace_strain_t]];

    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

    PylithScalar equivalentPressure_t = 0.0;

    // Equivalent pressure
    for (PylithInt i = 0; i < _phases; i++) {
        equivalentPressure_t += saturation[i]*pressure_t[i];
    }
    
    // Update porosity
    porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                       ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                       drainedBulkModulus * equivalentPressure_t);

} // updatePorosityImplicit

// ---------------------------------------------------------------------------------------------------------------------
/* Update saturation for a linear poroelastic material, implicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOilPlaneStrain::updateSaturationImplicit(const PylithInt dim,
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
                                                                                            PylithScalar saturation[]) {
    const PylithInt _dim = 2;
    const PylithInt _phases = 3;
    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_threePhaseSaturation = numA - 7;
    const PylithInt i_threePhaseFluidModulus = numA - 6;

    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA -2;

    // Constants
    // const PylithScalar dt = constants[0];

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosity] >= 0);
    assert(saturation);

    // Do stuff
    const PylithScalar* pressure = &s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];



    const PylithScalar* saturation_old = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_threePhaseFluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Generate Zeta
    PylithScalar Zeta[_phases], ZetaSum = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        Zeta[i] = 0.0;
    }

    // Nww
    N_tensor[0] = porosity*saturation_old[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[1];
    // Nwg
    N_tensor[2] = saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation_old[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation_old[1];
    // Noo
    N_tensor[4] = porosity*saturation_old[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation_old[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation_old[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2];
    // Ngw
    N_tensor[6] = saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[1];
    // Ngg
    N_tensor[8] = porosity*saturation_old[2]*(1.0 / fluidModulus[2]) + saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Zeta[i] += N_tensor[i*_phases+j]*pressure[j];
        } // for
        Zeta[i] += saturation_old[i]*biotCoefficient*trace_strain;
    } // for

    // Update saturation
    
    for (PylithInt i = 0; i < _phases; ++i ) {
        ZetaSum += Zeta[i];
    }
    for (PylithInt j = 0; j < _phases; j++) {
            saturation[j] += Zeta[j] / ZetaSum;
    }
} // updateSaturationImplicit


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
    const PylithInt _phases = 3;

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
    const PylithInt i_threePhaseSaturation = numA - 7;    
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

    // Do stuff
    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    PylithScalar equivalentPressure_t = 0.0;

    // Equivalent pressure
    for (PylithInt i = 0; i < _phases; i++) {
        equivalentPressure_t += saturation[i]*pressure_t[i];
    }

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * equivalentPressure_t);
} // updatePorosityExplicit


// =====================================================================================================================
// Kernels for isotropic, linear poroelasticity in 3D.
// =====================================================================================================================

// ----------------------------------------------------------------------

// ================================= STD =======================================

// ================================= LHS =======================================

// ---------------------------------------------------------------------------------------------------------------------
// f0u placeholder function for multiphaseporoelasticity equation
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0u_implicit(const PylithInt dim,
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
    // Incoming solution fields.

    // Incoming auxiliary fields.

    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += 0.0;
        // PetscPrintf(PETSC_COMM_WORLD, "f0u[%i]: %f\n",i, f0[i]);
    } // for
} // f0u

// ---------------------------------------------------------------------------------------------------------------------
// f0u placeholder function for multiphaseporoelasticity equation
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0u_body_implicit(const PylithInt dim,
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
    // Incoming solution fields.

    // Incoming auxiliary fields

    // MultiphasePoroelasticity

    // 3 + n
    const PylithInt i_bodyForce = 4;

    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += bodyForce[i];
    } // for
} // f0u_body_implicit

// ---------------------------------------------------------------------------------------------------------------------
// f0u_grav_implicit - f0 function for generic multiphaseporoelasticity terms ( + grav body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0u_grav_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;                                                          
    // Incoming auxililary fields.

    // MultiphasePoroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_porosity = 1;

    // 3 + n
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluid_density = numA - 13;
    const PylithInt i_threePhaseSaturation = numA - 7;

    const PylithScalar* fluidDensity = &a[aOff[i_fluid_density]];
    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar solidDensity = a[aOff[i_solid_density]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    PylithScalar equivalentFluidDensity = 0.0, bulkDensity = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        equivalentFluidDensity += saturation[i]*fluidDensity[i];
    }

    bulkDensity = (1.0 - porosity)*solidDensity + porosity*equivalentFluidDensity;

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += bulkDensity * gravityField[i];
    } // for

} // f0u_grav_implicit

// ---------------------------------------------------------------------------------------------------------------------
//f0u_body_grav_implicit - f0 function for generic multiphaseporoelasticity terms ( + grav body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0u_body_grav_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;                                                          
    // Incoming auxililary fields.

    // MultiphasePoroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_porosity = 1;

    // 3 + n
    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_fluid_density = numA - 13;
    const PylithInt i_threePhaseSaturation = numA - 7;

    const PylithScalar* fluidDensity = &a[aOff[i_fluid_density]];
    const PylithScalar* saturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar solidDensity = a[aOff[i_solid_density]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    PylithScalar equivalentFluidDensity = 0.0, bulkDensity = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        equivalentFluidDensity += saturation[i]*fluidDensity[i];
    }

    bulkDensity = (1.0 - porosity)*solidDensity + porosity*equivalentFluidDensity;

    // gravity field
    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += bulkDensity * gravityField[i];
    } // for

    // body force
    for (PylithInt i = 0; i < _dim; ++i) {
        f0[i] += bodyForce[i];
    } // for    

} // f0u_body_grav_implicit



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
    // const PylithInt _dim = 3;

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
    const PylithInt _phases = 3;

    // Incoming re-packed solution field.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t;
    } // for
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 4;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for

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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 5;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
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
    const PylithInt _phases = 3;

    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Poroelasticity
    const PylithInt i_porosity = 3;
    const PylithInt i_source = 6;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

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
    assert(f0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithScalar source = a[aOff[i_source]];

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            f0[i] += N_tensor[i*_phases+j]*pressure_t[j];
        } // for
        f0[i] += saturation[i]*biotCoefficient*trace_strain_t - saturation[i]*source;
    } // for
} // f0p_implicit_source_grav_body


// ---------------------------------------------------------------------------------------------------------------------
// f0pdot function for multiphaseporoelasticity equation, implicit time stepping, quasistatic.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::f0pdot_implicit(const PylithInt dim,
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
    const PylithInt _phases = 3;
    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_pdot = 4;

    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_pdot] >= 0);
    assert(s);
    assert(s_t);

    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]]; // pressure_t
    const PylithScalar* pdot = &s[sOff[i_pdot]]; // pdot

    for (PylithInt i = 0; i < _phases; ++i) {
        f0[i] += pressure_t[i];
        f0[i] -= pdot[i];
    }
} // f0pdot_implicit

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


    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
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
    //const PylithInt _phases = 3;

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_bodyForce = 4;

    // IsotropicLinearPoroelasticityBlackOil
    // const PylithInt i_fluidDensity = numA - 13;
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

    //const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - bodyForce[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - bodyForce[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[2] - bodyForce[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[3] - bodyForce[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[4] - bodyForce[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[5] - bodyForce[1]);

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
    //const PylithInt _phases = 3;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    //const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

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

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - bodyForce[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - bodyForce[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - bodyForce[j]);
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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;
    const PylithInt i_gravityField = 4;

    // IsotropicLinearPoroelasticityBlackOilBlackOil
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
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - fluidDensity[0] * gravityField[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - fluidDensity[0] * gravityField[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[0] - fluidDensity[1] * gravityField[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[1] - fluidDensity[1] * gravityField[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[0] - fluidDensity[2] * gravityField[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[1] - fluidDensity[2] * gravityField[1]);

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_gravityField = 4;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar *gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

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

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - fluidDensity[0] * gravityField[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - fluidDensity[1] * gravityField[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - fluidDensity[2] * gravityField[j]);
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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity
    const PylithInt i_fluidDensity = 1;
    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_threePhaseViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    vel_w[0] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[0] - bodyForce[0] - fluidDensity[0] * gravityField[0]);
    vel_w[1] = ( (isotropicPermeability * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[1] - bodyForce[1] - fluidDensity[0] * gravityField[1]);

    // Oil
    vel_o[0] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[0] - bodyForce[0] - fluidDensity[1] * gravityField[0]);
    vel_o[1] = ( (isotropicPermeability * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[1] - bodyForce[1] - fluidDensity[1] * gravityField[1]);

    // Gas
    vel_g[0] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[0] - bodyForce[0] - fluidDensity[2] * gravityField[0]);
    vel_g[1] = ( (isotropicPermeability * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[1] - bodyForce[1] - fluidDensity[2] * gravityField[1]);

    // Water
    f1[0] += vel_w[0] / formationVolumeFactor[0];
    f1[1] += vel_w[1] / formationVolumeFactor[0];

    // Oil
    f1[2] += vel_o[0] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[0]) / formationVolumeFactor[2];
    f1[3] += vel_o[1] / formationVolumeFactor[1] + (solutionOilGasRatio * vel_g[1]) / formationVolumeFactor[2];

    // Gas
    f1[4] += vel_g[0] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[0]) / formationVolumeFactor[1];
    f1[5] += vel_g[1] / formationVolumeFactor[2] + (solutionGasOilRatio * vel_o[1]) / formationVolumeFactor[1];

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
    //const PylithInt _phases = 3;

    // Incoming solution field.
    const PylithInt i_pressure = 1;

    // Incoming auxiliary field.

    // Poroelasticity

    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(f1);

    const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];

    const PylithScalar* fluidDensity = &a[aOff[i_fluidDensity]];
    const PylithScalar solutionOilGasRatio = a[aOff[i_solutionOilGasRatio]];
    const PylithScalar solutionGasOilRatio = a[aOff[i_solutionGasOilRatio]];
    const PylithScalar* formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar* relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar* fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

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

    PylithScalar vel_w[_dim], vel_o[_dim], vel_g[_dim];

    // Water
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_w[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[0]) / fluidViscosity[0]) * (pressure_x[j+0*_dim] - bodyForce[j] - fluidDensity[0] * gravityField[j]);
        } // for
    } // for

    // Oil
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_o[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[1]) / fluidViscosity[1]) * (pressure_x[j+1*_dim] - bodyForce[j] - fluidDensity[1] * gravityField[j]);
        } // for
    } // for

    // Gas
    for (PylithInt i = 0; i < _dim; ++i) {
        for (PylithInt j = 0; j < _dim; j++) {
            vel_g[i] += ( (tensorPermeability[i*_dim+j] * relativePermeability[2]) / fluidViscosity[2]) * (pressure_x[j+2*_dim] - bodyForce[j] - fluidDensity[2] * gravityField[j]);
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
    const PylithInt _phases = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_isotropicPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_isotropicPermeability] >= 0);
    assert(Jf3);

    const PylithScalar *formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar *relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar *fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];

    PylithScalar solutionRatios[_dim*_dim], coeff[_dim*_dim];
    PylithScalar tensorPermeability[_dim*_dim];

    tensorPermeability[0] = 0.0;
    tensorPermeability[1] = 0.0;
    tensorPermeability[2] = 0.0;
    tensorPermeability[3] = 0.0;
    tensorPermeability[4] = 0.0;
    tensorPermeability[5] = 0.0;
    tensorPermeability[6] = 0.0;
    tensorPermeability[7] = 0.0;
    tensorPermeability[8] = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        tensorPermeability[d*_dim+d] += isotropicPermeablity;
    }

    solutionRatios[0] = 1.0; // R_ww
    solutionRatios[1] = 0.0; // R_wo
    solutionRatios[2] = 0.0; // R_wg
    solutionRatios[3] = 0.0; // R_ow
    solutionRatios[4] = 1.0; // R_oo
    solutionRatios[5] = a[aOff[i_solutionOilGasRatio]]; // R_og
    solutionRatios[6] = 0.0; // R_gw
    solutionRatios[7] = a[aOff[i_solutionGasOilRatio]]; // R_go
    solutionRatios[8] = 1.0; // R_gg

    coeff[0] = (solutionRatios[0] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ww
    coeff[1] = (solutionRatios[1] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // wo
    coeff[2] = (solutionRatios[2] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // wg

    coeff[3] = (solutionRatios[3] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ow
    coeff[4] = (solutionRatios[4] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // oo
    coeff[5] = (solutionRatios[5] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // og

    coeff[6] = (solutionRatios[6] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // gw
    coeff[7] = (solutionRatios[7] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // go
    coeff[8] = (solutionRatios[8] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // gg

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            for (PylithInt k = 0; k < _dim; ++k) {
                for (PylithInt l = 0; l < _dim; ++l) {
                    Jf3[((i * _phases + j) * _dim + k) * _dim + l] -= coeff[i*_phases+j] * tensorPermeability[k*_dim+l];
                }
            }
        }
    }

    /* j(n_i,n_j,dim,dim) = \frac{\partial \vec{q}}{\partial \vec{p}}
     *
     *  0:  j0000 = ww_coeff * k00
     *  1:  j0001 = ww_coeff * k01
     *  2:  j0002 = ww_coeff * k02
     *  3:  j0010 = ww_coeff * k10
     *  4:  j0011 = ww_coeff * k11
     *  5:  j0012 = ww_coeff * k12
     *  6:  j0020 = ww_coeff * k20
     *  7:  j0021 = ww_coeff * k21
     *  8:  j0022 = ww_coeff * k22
     *  9:  j0100 = wo_coeff * k00
     * 10:  j0101 = wo_coeff * k01
     * 11:  j0102 = wo_coeff * k02
     * 12:  j0110 = wo_coeff * k10
     * 13:  j0111 = wo_coeff * k11
     * 14:  j0112 = wo_coeff * k12
     * 15:  j0120 = wo_coeff * k20
     * 16:  j0121 = wo_coeff * k21
     * 17:  j0122 = wo_coeff * k22
     * 18:  j0200 = wg_coeff * k00
     * 19:  j0201 = wg_coeff * k01
     * 20:  j0202 = wg_coeff * k02
     * 21:  j0210 = wg_coeff * k10
     * 22:  j0211 = wg_coeff * k11
     * 23:  j0212 = wg_coeff * k12
     * 24:  j0220 = wg_coeff * k20
     * 25:  j0221 = wg_coeff * k21
     * 26:  j0222 = wg_coeff * k22
     * 27:  j1000 = ow_coeff * k00
     * 28:  j1001 = ow_coeff * k01
     * 29:  j1002 = ow_coeff * k02
     * 30:  j1010 = ow_coeff * k10
     * 31:  j1011 = ow_coeff * k11
     * 32:  j1012 = ow_coeff * k12
     * 33:  j1020 = ow_coeff * k20
     * 34:  j1021 = ow_coeff * k21
     * 35:  j1022 = ow_coeff * k22
     * 36:  j1100 = oo_coeff * k00
     * 37:  j1101 = oo_coeff * k01
     * 38:  j1102 = oo_coeff * k02
     * 39:  j1110 = oo_coeff * k10
     * 40:  j1111 = oo_coeff * k11
     * 41:  j1112 = oo_coeff * k12
     * 42:  j1120 = oo_coeff * k20
     * 43:  j1121 = oo_coeff * k21
     * 44:  j1122 = oo_coeff * k22
     * 45:  j1200 = og_coeff * k00
     * 46:  j1201 = og_coeff * k01
     * 47:  j1202 = og_coeff * k02
     * 48:  j1210 = og_coeff * k10
     * 49:  j1211 = og_coeff * k11
     * 50:  j1212 = og_coeff * k12
     * 51:  j1220 = og_coeff * k20
     * 52:  j1221 = og_coeff * k21
     * 53:  j1222 = og_coeff * k22
     * 54:  j2000 = gw_coeff * k00
     * 55:  j2001 = gw_coeff * k01
     * 56:  j2002 = gw_coeff * k02
     * 57:  j2010 = gw_coeff * k10
     * 58:  j2011 = gw_coeff * k11
     * 59:  j2012 = gw_coeff * k12
     * 60:  j2020 = gw_coeff * k20
     * 61:  j2021 = gw_coeff * k21
     * 62:  j2022 = gw_coeff * k22
     * 63:  j2100 = go_coeff * k00
     * 64:  j2101 = go_coeff * k01
     * 65:  j2102 = go_coeff * k02
     * 66:  j2110 = go_coeff * k10
     * 67:  j2111 = go_coeff * k11
     * 68:  j2112 = go_coeff * k12
     * 69:  j2120 = go_coeff * k20
     * 70:  j2121 = go_coeff * k21
     * 71:  j2122 = go_coeff * k22
     * 72:  j2200 = gg_coeff * k00
     * 73:  j2201 = gg_coeff * k01
     * 74:  j2202 = gg_coeff * k02
     * 75:  j2210 = gg_coeff * k10
     * 76:  j2211 = gg_coeff * k11
     * 77:  j2212 = gg_coeff * k12
     * 78:  j2220 = gg_coeff * k20
     * 79:  j2221 = gg_coeff * k21
     * 80:  j2222 = gg_coeff * k22
     */

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
    const PylithInt _phases = 3;

    // index of Incoming auxiliary fields.
    // Poroelasticity
    const PylithInt i_fluidViscosity = 2;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_solutionOilGasRatio = numA - 12;
    const PylithInt i_solutionGasOilRatio = numA - 11;
    const PylithInt i_formationVolumeFactor = numA - 10;
    const PylithInt i_relativePermeability = numA - 9;
    const PylithInt i_threePhaseViscosity = numA - 8;
    const PylithInt i_tensorPermeability = numA - 1;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_fluidViscosity] >= 0);
    assert(aOff[i_tensorPermeability] >= 0);
    assert(Jf3);

    const PylithScalar *formationVolumeFactor = &a[aOff[i_formationVolumeFactor]];
    const PylithScalar *relativePermeability = &a[aOff[i_relativePermeability]];
    const PylithScalar *fluidViscosity = &a[aOff[i_threePhaseViscosity]];
    const PylithScalar* vectorPermeability = &a[aOff[i_tensorPermeability]];

    PylithScalar solutionRatios[_dim*_dim], coeff[_dim*_dim];
    PylithScalar tensorPermeability[_dim*_dim];

    tensorPermeability[0] = vectorPermeability[0];
    tensorPermeability[1] = vectorPermeability[3];
    tensorPermeability[2] = vectorPermeability[5];
    tensorPermeability[3] = vectorPermeability[3];
    tensorPermeability[4] = vectorPermeability[1];
    tensorPermeability[5] = vectorPermeability[4];
    tensorPermeability[6] = vectorPermeability[5];
    tensorPermeability[7] = vectorPermeability[4];
    tensorPermeability[8] = vectorPermeability[2];

    solutionRatios[0] = 1.0; // R_ww
    solutionRatios[1] = 0.0; // R_wo
    solutionRatios[2] = 0.0; // R_wg
    solutionRatios[3] = 0.0; // R_ow
    solutionRatios[4] = 1.0; // R_oo
    solutionRatios[5] = a[aOff[i_solutionOilGasRatio]]; // R_og
    solutionRatios[6] = 0.0; // R_gw
    solutionRatios[7] = a[aOff[i_solutionGasOilRatio]]; // R_go
    solutionRatios[8] = 1.0; // R_gg

    coeff[0] = (solutionRatios[0] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ww
    coeff[1] = (solutionRatios[1] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // wo
    coeff[2] = (solutionRatios[2] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // wg

    coeff[3] = (solutionRatios[3] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // ow
    coeff[4] = (solutionRatios[4] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // oo
    coeff[5] = (solutionRatios[5] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // og

    coeff[6] = (solutionRatios[6] * relativePermeability[0]) / (formationVolumeFactor[0]*fluidViscosity[0]); // gw
    coeff[7] = (solutionRatios[7] * relativePermeability[1]) / (formationVolumeFactor[1]*fluidViscosity[1]); // go
    coeff[8] = (solutionRatios[8] * relativePermeability[2]) / (formationVolumeFactor[2]*fluidViscosity[2]); // gg

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            for (PylithInt k = 0; k < _dim; ++k) {
                for (PylithInt l = 0; l < _dim; ++l) {
                    Jf3[((i * _phases + j) * _dim + k) * _dim + l] -= coeff[i*_phases+j] * tensorPermeability[k*_dim+l];
                }
            }
        }
    }

    /* j(n_i,n_j,dim,dim) = \frac{\partial \vec{q}}{\partial \vec{p}}
     *
     *  0:  j0000 = ww_coeff * k00
     *  1:  j0001 = ww_coeff * k01
     *  2:  j0002 = ww_coeff * k02
     *  3:  j0010 = ww_coeff * k10
     *  4:  j0011 = ww_coeff * k11
     *  5:  j0012 = ww_coeff * k12
     *  6:  j0020 = ww_coeff * k20
     *  7:  j0021 = ww_coeff * k21
     *  8:  j0022 = ww_coeff * k22
     *  9:  j0100 = wo_coeff * k00
     * 10:  j0101 = wo_coeff * k01
     * 11:  j0102 = wo_coeff * k02
     * 12:  j0110 = wo_coeff * k10
     * 13:  j0111 = wo_coeff * k11
     * 14:  j0112 = wo_coeff * k12
     * 15:  j0120 = wo_coeff * k20
     * 16:  j0121 = wo_coeff * k21
     * 17:  j0122 = wo_coeff * k22
     * 18:  j0200 = wg_coeff * k00
     * 19:  j0201 = wg_coeff * k01
     * 20:  j0202 = wg_coeff * k02
     * 21:  j0210 = wg_coeff * k10
     * 22:  j0211 = wg_coeff * k11
     * 23:  j0212 = wg_coeff * k12
     * 24:  j0220 = wg_coeff * k20
     * 25:  j0221 = wg_coeff * k21
     * 26:  j0222 = wg_coeff * k22
     * 27:  j1000 = ow_coeff * k00
     * 28:  j1001 = ow_coeff * k01
     * 29:  j1002 = ow_coeff * k02
     * 30:  j1010 = ow_coeff * k10
     * 31:  j1011 = ow_coeff * k11
     * 32:  j1012 = ow_coeff * k12
     * 33:  j1020 = ow_coeff * k20
     * 34:  j1021 = ow_coeff * k21
     * 35:  j1022 = ow_coeff * k22
     * 36:  j1100 = oo_coeff * k00
     * 37:  j1101 = oo_coeff * k01
     * 38:  j1102 = oo_coeff * k02
     * 39:  j1110 = oo_coeff * k10
     * 40:  j1111 = oo_coeff * k11
     * 41:  j1112 = oo_coeff * k12
     * 42:  j1120 = oo_coeff * k20
     * 43:  j1121 = oo_coeff * k21
     * 44:  j1122 = oo_coeff * k22
     * 45:  j1200 = og_coeff * k00
     * 46:  j1201 = og_coeff * k01
     * 47:  j1202 = og_coeff * k02
     * 48:  j1210 = og_coeff * k10
     * 49:  j1211 = og_coeff * k11
     * 50:  j1212 = og_coeff * k12
     * 51:  j1220 = og_coeff * k20
     * 52:  j1221 = og_coeff * k21
     * 53:  j1222 = og_coeff * k22
     * 54:  j2000 = gw_coeff * k00
     * 55:  j2001 = gw_coeff * k01
     * 56:  j2002 = gw_coeff * k02
     * 57:  j2010 = gw_coeff * k10
     * 58:  j2011 = gw_coeff * k11
     * 59:  j2012 = gw_coeff * k12
     * 60:  j2020 = gw_coeff * k20
     * 61:  j2021 = gw_coeff * k21
     * 62:  j2022 = gw_coeff * k22
     * 63:  j2100 = go_coeff * k00
     * 64:  j2101 = go_coeff * k01
     * 65:  j2102 = go_coeff * k02
     * 66:  j2110 = go_coeff * k10
     * 67:  j2111 = go_coeff * k11
     * 68:  j2112 = go_coeff * k12
     * 69:  j2120 = go_coeff * k20
     * 70:  j2121 = go_coeff * k21
     * 71:  j2122 = go_coeff * k22
     * 72:  j2200 = gg_coeff * k00
     * 73:  j2201 = gg_coeff * k01
     * 74:  j2202 = gg_coeff * k02
     * 75:  j2210 = gg_coeff * k10
     * 76:  j2211 = gg_coeff * k11
     * 77:  j2212 = gg_coeff * k12
     * 78:  j2220 = gg_coeff * k20
     * 79:  j2221 = gg_coeff * k21
     * 80:  j2222 = gg_coeff * k22
     */

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
    const PylithInt _phases = 3;

    // Incoming auxiliary fields.

        // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(Jf0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Jf0[i*_phases + j] = utshift * N_tensor[i*_phases + j];
        }
    }
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
    const PylithInt _phases = 3;

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

    for (PylithInt i = 0; i < _phases; ++i) { 
    Jf0[i] += utshift * biotCoefficient;
    }
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
    const PylithInt _phases = 3;

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_saturation = numA - 7;
    const PylithInt i_fluidModulus = numA - 6;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA - 2;

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 2);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_biotCoefficient] >= 0);
    assert(Jf0);

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithScalar* saturation = &a[aOff[i_saturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_fluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Nww
    N_tensor[0] = porosity*saturation[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Nwg
    N_tensor[2] = saturation[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation[1];
    // Noo
    N_tensor[4] = porosity*saturation[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2];
    // Ngw
    N_tensor[6] = saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[1];
    // Ngg
    N_tensor[8] = porosity*saturation[2]*(1.0 / fluidModulus[2]) + saturation[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Jf0[i*_phases + j] = N_tensor[i*_phases + j];
        }
    }
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
    const PylithInt _phases = 3;

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

    for (PylithInt i = 0; i < _phases; ++i) {
        Jf0[i] += biotCoefficient;
    }
} // Jf0pedot

// -----------------------------------------------------------------------------
// Jf0pdotp - Jf0 function for isotropic linear multiphaseporoelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0pdotp(const PylithInt dim,
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
                                                      PylithScalar Jf0[]) {
    const PylithInt _phases = 3;
    assert(aOff);
    assert(a);

    for (PetscInt i = 0; i < _phases; ++i) {
        Jf0[i*_phases+i] += s_tshift;
    }
} // Jf0pdotp


// -----------------------------------------------------------------------------
// Jf0pdotpdot - Jf0 function for isotropic linear multiphaseporoelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::Jf0pdotpdot(const PylithInt dim,
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
                                                         PylithScalar Jf0[]) {
    const PylithInt _phases = 3;
    assert(aOff);
    assert(a);

    for (PetscInt i = 0; i < _phases; ++i) {
        Jf0[i*_phases+i] -= 1.0;
    }
} // Jf0pdotpdot


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
    // const PylithInt i_pressure = 1;
    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 4;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 5;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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

    const PylithInt i_velocity = 2;

    // Incoming re-packed auxiliary field.
    // Poroelasticity
    const PylithInt i_source = 6;
    const PylithScalar source = a[aOff[i_source]];

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];


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
    //const PylithInt i_velocity = 2;

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

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    //const PylithScalar biotModulus = a[aOff[i_biotModulus]];

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
    const PylithInt _phases = 3;

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
    const PylithInt i_threePhaseSaturation = numA - 7;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

    // Do stuff
    const PylithScalar* pressure_t = &s[sOff[i_pressure_t]];
    const PylithScalar trace_strain_t = s[sOff[i_trace_strain_t]];

    const PylithScalar* threePhaseSaturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar equivalentPressure_t = 0.0;

    // Equivalent pressure
    for (PylithInt i = 0; i < _phases; i++) {
        equivalentPressure_t += threePhaseSaturation[i]*pressure_t[i];
    }

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * equivalentPressure_t);
} // updatePorosityImplicit

// ---------------------------------------------------------------------------------------------------------------------
/* Update saturation for a linear poroelastic material, implicit.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityBlackOil3D::updateSaturationImplicit(const PylithInt dim,
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
                                                                                            PylithScalar saturation[]) {
    const PylithInt _dim = 3;
    const PylithInt _phases = 3;
    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_trace_strain = 2;

    // Incoming re-packed auxiliary field.

    // Poroelasticity
    const PylithInt i_porosity = 3;

    // IsotropicLinearPoroelasticityBlackOil
    const PylithInt i_threePhaseSaturation = numA - 7;
    const PylithInt i_threePhaseFluidModulus = numA - 6;

    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_solidBulkModulus = numA -2;

    // Constants
    //const PylithScalar dt = constants[0];

    // Run Checks
    assert(_dim == dim);
    assert(numS >= 3);
    assert(numA >= 3);
    assert(aOff);
    assert(aOff[i_porosity] >= 0);
    assert(saturation);

    // Do stuff
    const PylithScalar* pressure = &s[sOff[i_pressure]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    //const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];

    const PylithScalar* saturation_old = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar* fluidModulus = &a[aOff[i_threePhaseFluidModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar solidBulkModulus = a[aOff[i_solidBulkModulus]];
    const PylithScalar porosity = a[aOff[i_porosity]];

    // Account for capiliarity
    PylithScalar dSw_dpco = 0.0;
    PylithScalar dSg_dpcg = 0.0;

    // Generate N Tensor
    PylithScalar N_tensor[_phases*_phases];

    // Generate Zeta
    PylithScalar Zeta[_phases], ZetaSum = 0.0;

    for (PylithInt i = 0; i < _phases; ++i) {
        Zeta[i] = 0.0;
    }

    // Nww
    N_tensor[0] = porosity*saturation_old[0] * (1.0 / fluidModulus[0]) - porosity*dSw_dpco + saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Nwo
    N_tensor[1] = porosity*dSw_dpco + saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[1];
    // Nwg
    N_tensor[2] = saturation_old[0]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2];
    // Now
    N_tensor[3] = porosity*dSw_dpco + saturation_old[1] * ( (biotCoefficient - porosity) / solidBulkModulus ) * saturation_old[1];
    // Noo
    N_tensor[4] = porosity*saturation_old[1] * (1.0 / fluidModulus[1]) + porosity*( -dSw_dpco + dSg_dpcg) + saturation_old[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Nog
    N_tensor[5] = porosity*(-dSg_dpcg) + saturation_old[1]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2];
    // Ngw
    N_tensor[6] = saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[0];
    // Ngo
    N_tensor[7] = -porosity*dSg_dpcg + saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[1];
    // Ngg
    N_tensor[8] = porosity*saturation_old[2]*(1.0 / fluidModulus[2]) + saturation_old[2]*( (biotCoefficient - porosity) / solidBulkModulus )*saturation_old[2] + porosity*dSg_dpcg;

    for (PylithInt i = 0; i < _phases; ++i) {
        for (PylithInt j = 0; j < _phases; ++j) {
            Zeta[i] += N_tensor[i*_phases+j]*pressure[j];
        } // for
        Zeta[i] += saturation_old[i]*biotCoefficient*trace_strain;
    } // for

    // Update saturation
    
    for (PylithInt i = 0; i < _phases; ++i ) {
        ZetaSum += Zeta[i];
    }
    for (PylithInt j = 0; j < _phases; j++) {
            saturation[j] += Zeta[j] / ZetaSum;
    }
} // updateSaturationImplicit


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
    const PylithInt _phases = 3;

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
    const PylithInt i_threePhaseSaturation = numA - 7;
    const PylithInt i_drainedBulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    // Constants
    const PylithScalar dt = constants[0];

    // Do stuff
    const PylithScalar* pressure_t = &s_t[sOff[i_pressure]];
    const PylithScalar* velocity_x = &s_x[sOff[i_velocity]];

    const PylithScalar* threePhaseSaturation = &a[aOff[i_threePhaseSaturation]];
    const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithScalar trace_strain_t = 0.0;

    for (PylithInt d = 0; d < _dim; ++d) {
        trace_strain_t += velocity_x[d*_dim+d];
    }

    PylithScalar equivalentPressure_t = 0.0;

    // Equivalent pressure
    for (PylithInt i = 0; i < _phases; i++) {
        equivalentPressure_t += threePhaseSaturation[i]*pressure_t[i];
    }

    // Update porosity
    porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                              ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                              drainedBulkModulus * equivalentPressure_t);
} // updatePorosityExplicit


// End of file
