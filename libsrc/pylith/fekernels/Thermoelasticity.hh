/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

/*
 * Kernels for thermoelasticity.
 *
 * Thermoelasticity couples:
 * 1. Momentum equation: ∇·σ + f = ρ∂²u/∂t² (or = 0 for quasistatic)
 * 2. Heat equation: ρc ∂T/∂t - ∇·(k∇T) = Q
 *
 * With constitutive relation including thermal strain:
 * σ = C:(ε - ε_thermal) = C:(ε - α(T-T_ref)I)
 *
 * Solution fields: [displacement, temperature]
 *
 * Auxiliary fields:
 *   Index 0: density
 *   Index 1: body_force (optional)
 *   Index 2: gravitational_acceleration (optional)
 *   Index 3: specific_heat
 *   Index 4: thermal_conductivity
 *   Index 5: heat_source (optional)
 *   Index 6: reference_temperature
 *   Index 7: thermal_expansion_coefficient
 *   Index 8+: rheology-specific (shear_modulus, bulk_modulus, etc.)
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert>

class pylith::fekernels::Thermoelasticity {
public:

    // ============================= Solution Field Indices =============================
    static const PylithInt i_displacement = 0;
    static const PylithInt i_temperature = 1;

    // ============================= Auxiliary Field Indices =============================
    static const PylithInt i_density = 0;
    static const PylithInt i_specificHeat = 1;
    static const PylithInt i_thermalConductivity = 2;
    // Optional fields shift indices for rheology fields
    static const PylithInt i_refTemperature = 3;
    static const PylithInt i_thermalExpansionCoeff = 4;
    // Rheology fields start at index 5 for the base case

    // ============================= Displacement Equation: LHS Residual ================================

    // ----------------------------------------------------------------------
    /** f0 function for displacement equation with inertia (LHS, dynamic).
     *
     * f0_u = ρ ∂²u/∂t² = ρ u_tt
     *
     * Solution fields: [displacement, temperature]
     * Auxiliary fields: [density, ...]
     */
    static inline
    void f0u_inertia(const PylithInt dim,
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
        assert(numS >= 2);
        assert(numA >= 1);
        assert(sOff[i_displacement] >= 0);
        assert(aOff[i_density] >= 0);
        assert(s_t);

        const PylithScalar* disp_tt = &s_t[sOff[i_displacement]]; // Using velocity from s_t
        const PylithScalar density = a[aOff[i_density]];

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += density * disp_tt[i];
        } // for
    } // f0u_inertia

    // ============================= Temperature Equation: LHS Residual ================================

    // ----------------------------------------------------------------------
    /** f0 function for temperature equation, time-dependent (LHS).
     *
     * f0_T = ρc ∂T/∂t
     *
     * Solution fields: [displacement, temperature]
     * Auxiliary fields: [density, specific_heat, ...]
     */
    static inline
    void f0T_timedep(const PylithInt dim,
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
        assert(numS >= 2);
        assert(numA >= 2);
        assert(sOff[i_temperature] >= 0);
        assert(aOff[i_density] >= 0);
        assert(aOff[i_specificHeat] >= 0);
        assert(s_t);

        const PylithScalar temperature_t = s_t[sOff[i_temperature]];
        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        f0[0] += density * specificHeat * temperature_t;
    } // f0T_timedep


    // ============================= Two-Way Thermoelastic Coupling ================================
    // For dynamic problems with solution [displacement, velocity, temperature]

    // Solution field indices for dynamic formulation
    static const PylithInt i_velocity = 1;
    static const PylithInt i_temperature_dynamic = 2;

    // ----------------------------------------------------------------------
    /** f0 function for temperature equation with thermoelastic heating.
     *
     * Two-way coupling: f0_T = ρc ∂T/∂t - T α 3K ∇·v
     *
     * The thermoelastic heating term T α 3K ∇·v represents adiabatic 
     * heating/cooling due to volumetric strain rate.
     * - Compression (∇·v < 0): heating
     * - Expansion (∇·v > 0): cooling
     *
     * Solution fields: [displacement, velocity, temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, 
     *                    reference_temperature, thermal_expansion_coefficient,
     *                    shear_modulus, bulk_modulus]
     */
    static inline
    void f0T_thermoelastic_twoway(const PylithInt dim,
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
        assert(numS >= 3);
        assert(numA >= 7);
        assert(s_t);

        // Auxiliary field indices
        const PylithInt i_thermalExpCoeff = 4;
        const PylithInt i_bulkMod = 6;

        const PylithScalar temperature = s[sOff[i_temperature_dynamic]];
        const PylithScalar temperature_t = s_t[sOff[i_temperature_dynamic]];
        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpCoeff]];
        const PylithScalar bulkModulus = a[aOff[i_bulkMod]];

        // Compute volumetric strain rate: ∇·v = ∂v_i/∂x_i
        const PylithScalar* vel_x = &s_x[sOff_x[i_velocity]];
        PylithScalar div_velocity = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            div_velocity += vel_x[d * dim + d];
        }

        // Heat capacity term: ρc ∂T/∂t
        f0[0] += density * specificHeat * temperature_t;

        // Thermoelastic heating term: -T α 3K ∇·v
        // Negative because heating occurs on compression (∇·v < 0)
        f0[0] -= temperature * thermalExpansionCoeff * 3.0 * bulkModulus * div_velocity;
    } // f0T_thermoelastic_twoway

    // ----------------------------------------------------------------------
    /** g0 function for thermoelastic heating source (RHS, explicit).
     *
     * For explicit time stepping, thermoelastic heating on RHS:
     * g0_T = T α 3K ∇·v
     *
     * Solution fields: [displacement, velocity, temperature]
     */
    static inline
    void g0T_thermoelastic(const PylithInt dim,
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
        assert(numS >= 3);
        assert(numA >= 7);

        const PylithInt i_thermalExpCoeff = 4;
        const PylithInt i_bulkMod = 6;

        const PylithScalar temperature = s[sOff[i_temperature_dynamic]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpCoeff]];
        const PylithScalar bulkModulus = a[aOff[i_bulkMod]];

        // Compute ∇·v
        const PylithScalar* vel_x = &s_x[sOff_x[i_velocity]];
        PylithScalar div_velocity = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            div_velocity += vel_x[d * dim + d];
        }

        // Thermoelastic heating: T α 3K ∇·v
        g0[0] += temperature * thermalExpansionCoeff * 3.0 * bulkModulus * div_velocity;
    } // g0T_thermoelastic

    // ============================= LHS Jacobian ================================

    // ----------------------------------------------------------------------
    /** Jf0_uu function for displacement equation with inertia.
     *
     * Jf0_uu = ρ * s_tshift
     */
    static inline
    void Jf0uu_inertia(const PylithInt dim,
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
        assert(numA >= 1);
        assert(aOff[i_density] >= 0);

        const PylithScalar density = a[aOff[i_density]];

        for (PylithInt i = 0; i < dim; ++i) {
            Jf0[i * dim + i] += density * s_tshift;
        } // for
    } // Jf0uu_inertia

    // ----------------------------------------------------------------------
    /** Jf0_TT function for temperature equation (time derivative).
     *
     * Jf0_TT = ρc * s_tshift
     */
    static inline
    void Jf0TT_timedep(const PylithInt dim,
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
        assert(numA >= 2);
        assert(aOff[i_density] >= 0);
        assert(aOff[i_specificHeat] >= 0);

        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        Jf0[0] += density * specificHeat * s_tshift;
    } // Jf0TT_timedep

    // ============================= Two-Way Coupling Jacobians ================================

    // ----------------------------------------------------------------------
    /** Jf0_TT for thermoelastic coupling (includes thermoelastic effect).
     *
     * Jf0_TT = ρc * s_tshift - α 3K ∇·v
     *
     * The second term comes from ∂(Tα3K∇·v)/∂T = α 3K ∇·v
     */
    static inline
    void Jf0TT_thermoelastic(const PylithInt dim,
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
        assert(numS >= 3);
        assert(numA >= 7);

        const PylithInt i_thermalExpCoeff = 4;
        const PylithInt i_bulkMod = 6;

        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpCoeff]];
        const PylithScalar bulkModulus = a[aOff[i_bulkMod]];

        // Compute ∇·v
        const PylithScalar* vel_x = &s_x[sOff_x[i_velocity]];
        PylithScalar div_velocity = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            div_velocity += vel_x[d * dim + d];
        }

        // Heat capacity contribution
        Jf0[0] += density * specificHeat * s_tshift;

        // Thermoelastic contribution: -α 3K ∇·v (derivative of -Tα3K∇·v w.r.t. T)
        Jf0[0] -= thermalExpansionCoeff * 3.0 * bulkModulus * div_velocity;
    } // Jf0TT_thermoelastic

    // ----------------------------------------------------------------------
    /** Jf3_Tv for thermoelastic coupling.
     *
     * Jf3_Tv = -T α 3K δ_ij (derivative of thermoelastic term w.r.t. velocity gradient)
     *
     * The thermoelastic heating is f0_T = ... - T α 3K ∇·v = ... - T α 3K (∂v_i/∂x_i)
     * So ∂f0_T/∂(∂v_j/∂x_k) = -T α 3K δ_jk
     *
     * This couples the temperature equation to velocity.
     */
    static inline
    void Jf3Tv_thermoelastic(const PylithInt dim,
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
        assert(numS >= 3);
        assert(numA >= 7);

        const PylithInt i_thermalExpCoeff = 4;
        const PylithInt i_bulkMod = 6;

        const PylithScalar temperature = s[sOff[i_temperature_dynamic]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpCoeff]];
        const PylithScalar bulkModulus = a[aOff[i_bulkMod]];

        const PylithScalar coupling = -temperature * thermalExpansionCoeff * 3.0 * bulkModulus;

        // Jf3[j*dim + k] = ∂f0_T/∂(∂v_j/∂x_k) = -T α 3K δ_jk
        for (PylithInt d = 0; d < dim; ++d) {
            Jf3[d * dim + d] += coupling;
        }
    } // Jf3Tv_thermoelastic

}; // class Thermoelasticity

// End of file
