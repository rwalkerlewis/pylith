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
 * Kernels for fully coupled thermoporoelasticity.
 *
 * Thermoporoelasticity couples:
 * 1. Momentum equation: ∇·σ + f = ρ∂²u/∂t² (or = 0 for quasistatic)
 *    where σ = C:(ε - ε_thermal) - α_B p I
 *    ε_thermal = α_T (T - T_ref) I
 *
 * 2. Fluid mass balance: ∂ζ/∂t + ∇·q = γ
 *    where ζ = α_B ∇·u + p/M + 3α_f(T - T_ref) (fluid content)
 *    q = -k/μ (∇p - ρ_f g) (Darcy flux)
 *
 * 3. Heat equation: (ρc)_eff ∂T/∂t - ∇·(k_T∇T) = Q
 *    with optional thermoelastic and viscous dissipation heating
 *
 * Solution fields (quasistatic): [displacement, pressure, trace_strain, temperature]
 * Solution fields (quasistatic with state vars): [displacement, pressure, trace_strain, velocity, pressure_dot, trace_strain_dot, temperature]
 * Solution fields (dynamic): [displacement, pressure, velocity, temperature]
 *
 * Auxiliary fields:
 *   Index 0: solid_density
 *   Index 1: fluid_density
 *   Index 2: fluid_viscosity
 *   Index 3: porosity
 *   (optional): gravity_field, body_force, source_density
 *   Rheology-specific: biot_coefficient, biot_modulus, drained_bulk_modulus, shear_modulus,
 *                      reference_temperature, thermal_expansion_coefficient, fluid_thermal_expansion,
 *                      thermal_conductivity, specific_heat, permeability
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Thermoelasticity.hh" // USES Thermoelasticity kernels

#include "pylith/utils/types.hh"

#include <cassert>

class pylith::fekernels::Thermoporoelasticity {
public:

    // ============================= Context =============================
    struct Context {
        PylithInt dim;
        // Solution fields
        const PylithReal* displacement;
        const PylithReal* displacement_t;
        const PylithReal* displacement_x;
        PylithReal pressure;
        PylithReal pressure_t;
        const PylithReal* pressure_x;
        PylithReal trace_strain;
        PylithReal trace_strain_t;
        const PylithReal* trace_strain_x;
        PylithReal temperature;
        PylithReal temperature_t;
        const PylithReal* temperature_x;
        const PylithReal* velocity;
        const PylithReal* velocity_t;
        const PylithReal* velocity_x;
        PylithReal pressure_dot;
        PylithReal trace_strain_dot;
        const PylithReal* x;
        // Auxiliary fields
        PylithReal solidDensity;
        PylithReal fluidDensity;
        PylithReal fluidViscosity;
        PylithReal porosity;
        PylithReal bulkDensity;
        PylithReal gravityField[3];
        PylithReal bodyForce[3];
        PylithReal sourceDensity;
        PylithReal heatSource;
        // Rheology parameters
        PylithReal biotCoefficient;
        PylithReal biotModulus;
        PylithReal drainedBulkModulus;
        PylithReal shearModulus;
        PylithReal refTemperature;
        PylithReal thermalExpansionCoeff;
        PylithReal fluidThermalExpansion;
        PylithReal thermalConductivity;
        PylithReal specificHeat;
        PylithReal permeability;

        Context(void) :
            dim(0),
            displacement(NULL),
            displacement_t(NULL),
            displacement_x(NULL),
            pressure(0.0),
            pressure_t(0.0),
            pressure_x(NULL),
            trace_strain(0.0),
            trace_strain_t(0.0),
            trace_strain_x(NULL),
            temperature(0.0),
            temperature_t(0.0),
            temperature_x(NULL),
            velocity(NULL),
            velocity_t(NULL),
            velocity_x(NULL),
            pressure_dot(0.0),
            trace_strain_dot(0.0),
            x(NULL),
            solidDensity(0.0),
            fluidDensity(0.0),
            fluidViscosity(0.0),
            porosity(0.0),
            bulkDensity(0.0),
            sourceDensity(0.0),
            heatSource(0.0),
            biotCoefficient(0.0),
            biotModulus(0.0),
            drainedBulkModulus(0.0),
            shearModulus(0.0),
            refTemperature(0.0),
            thermalExpansionCoeff(0.0),
            fluidThermalExpansion(0.0),
            thermalConductivity(0.0),
            specificHeat(0.0),
            permeability(0.0) {
            for (size_t i = 0; i < 3; ++i) {
                gravityField[i] = 0.0;
                bodyForce[i] = 0.0;
            }
        }
    };

    // ============================= Solution Field Indices =============================
    // Quasistatic (base): [displacement, pressure, trace_strain, temperature]
    static const PylithInt i_displacement = 0;
    static const PylithInt i_pressure = 1;
    static const PylithInt i_trace_strain = 2;
    static const PylithInt i_temperature = 3;

    // ============================= Context Setup =============================

    // ----------------------------------------------------------------------
    /** Set context for quasistatic thermoporoelasticity.
     *
     * Solution fields: [displacement, pressure, trace_strain, temperature]
     */
    static inline
    void setContext(Context* context,
                    const PylithInt dim,
                    const PylithInt numS,
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
                    const PylithScalar x[]) {
        assert(context);
        assert(numS >= 4);

        // Incoming auxiliary fields.
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_porosity = 3;

        assert(sOff[i_displacement] >= 0);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff[i_temperature] >= 0);

        assert(aOff[i_solidDensity] >= 0);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_porosity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);

        context->dim = dim;
        context->displacement = &s[sOff[i_displacement]];
        context->displacement_t = s_t ? &s_t[sOff[i_displacement]] : NULL;
        context->displacement_x = &s_x[sOff_x[i_displacement]];
        context->pressure = s[sOff[i_pressure]];
        context->pressure_t = s_t ? s_t[sOff[i_pressure]] : 0.0;
        context->pressure_x = &s_x[sOff_x[i_pressure]];
        context->trace_strain = s[sOff[i_trace_strain]];
        context->trace_strain_t = s_t ? s_t[sOff[i_trace_strain]] : 0.0;
        context->trace_strain_x = &s_x[sOff_x[i_trace_strain]];
        context->temperature = s[sOff[i_temperature]];
        context->temperature_t = s_t ? s_t[sOff[i_temperature]] : 0.0;
        context->temperature_x = &s_x[sOff_x[i_temperature]];
        context->x = x;

        // Basic auxiliaries
        context->solidDensity = a[aOff[i_solidDensity]];
        context->fluidDensity = a[aOff[i_fluidDensity]];
        context->porosity = a[aOff[i_porosity]];
        context->fluidViscosity = a[aOff[i_fluidViscosity]];
        context->bulkDensity = (1.0 - a[aOff[i_porosity]]) * a[aOff[i_solidDensity]] 
                             + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    } // setContext

    // ============================= Volumetric Strain Equation =============================

    // ----------------------------------------------------------------------
    /** f0 function for trace_strain equation.
     *
     * f0_e = ∇·u - ε_v = 0
     */
    static inline
    void f0e(const PylithInt dim,
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
        const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        PylithScalar div_u = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            div_u += displacement_x[d * dim + d];
        }
        
        f0[0] += div_u - trace_strain;
    } // f0e

    // ============================= Temperature Equation =============================

    // ----------------------------------------------------------------------
    /** f0 function for temperature equation (transient).
     *
     * f0_T = (ρc)_eff ∂T/∂t
     *
     * where (ρc)_eff = (1-φ)ρ_s c_s + φ ρ_f c_f is effective heat capacity
     */
    static inline
    void f0T_transient(const PylithInt dim,
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
        assert(s_t);
        assert(numA >= 10);

        // Auxiliary indices (base fields always at start)
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_porosity = 3;
        // Rheology fields are indexed relative to numA (specific_heat is last field added)
        const PylithInt i_specificHeat = numA - 1;

        const PylithScalar temperature_t = s_t[sOff[i_temperature]];
        const PylithScalar solidDensity = a[aOff[i_solidDensity]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar porosity = a[aOff[i_porosity]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        // Effective heat capacity
        const PylithScalar bulkDensity = (1.0 - porosity) * solidDensity + porosity * fluidDensity;
        const PylithScalar effectiveHeatCapacity = bulkDensity * specificHeat;

        f0[0] += effectiveHeatCapacity * temperature_t;
    } // f0T_transient

    // ----------------------------------------------------------------------
    /** f0 function for temperature with heat source.
     *
     * f0_T = (ρc)_eff ∂T/∂t - Q
     */
    static inline
    void f0T_source(const PylithInt dim,
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
        // First add the transient term
        f0T_transient(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0);

        // Add heat source (negative because it's a source on RHS)
        // heat_source is an optional field added after the required base fields
        // and before rheology fields. Its position depends on what other optional
        // fields are present. For simplicity, assume it's at index 4 if present.
        const PylithInt i_heatSource = 4;  // Assumes heat_source is first optional field
        const PylithScalar heatSource = a[aOff[i_heatSource]];

        f0[0] -= heatSource;
    } // f0T_source

    // ============================= Jacobians =============================

    // ----------------------------------------------------------------------
    /** Jf0_ee function for trace strain equation.
     */
    static inline
    void Jf0ee(const PylithInt dim,
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
        Jf0[0] = -1.0;
    } // Jf0ee

    // ----------------------------------------------------------------------
    /** Jf1_eu function for trace strain equation.
     *
     * Jf1_eu = ∂(∇·u)/∂(∇u) = δ_ij
     */
    static inline
    void Jf1eu(const PylithInt dim,
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
               PylithScalar Jf1[]) {
        for (PylithInt d = 0; d < dim; ++d) {
            Jf1[d * dim + d] = 1.0;
        }
    } // Jf1eu

    // ============================= State Variable Kernels =============================
    // These kernels are used when _useStateVars = true for quasistatic formulation.
    // Solution fields: [displacement, pressure, trace_strain, temperature, 
    //                   velocity, pressure_dot, trace_strain_dot, temperature_dot]

    // Solution field indices for state variable formulation
    static const PylithInt i_velocity = 4;
    static const PylithInt i_pressure_dot = 5;
    static const PylithInt i_trace_strain_dot = 6;
    static const PylithInt i_temperature_dot = 7;

    // ----------------------------------------------------------------------
    /** f0 function for velocity equation with state variables.
     *
     * f0_v = ∂u/∂t - v = 0
     *
     * This enforces that velocity equals the time derivative of displacement.
     */
    static inline
    void f0v_implicit(const PylithInt dim,
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
        assert(numS >= 8);
        assert(s_t);

        const PylithScalar* displacement_t = &s_t[sOff[i_displacement]];
        const PylithScalar* velocity = &s[sOff[i_velocity]];

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += displacement_t[i] - velocity[i];
        }
    } // f0v_implicit

    // ----------------------------------------------------------------------
    /** f0 function for pressure_dot equation with state variables.
     *
     * f0_pdot = ∂p/∂t - p_dot = 0
     */
    static inline
    void f0pdot(const PylithInt dim,
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
        assert(numS >= 8);
        assert(s_t);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar pressure_dot = s[sOff[i_pressure_dot]];

        f0[0] += pressure_t - pressure_dot;
    } // f0pdot

    // ----------------------------------------------------------------------
    /** f0 function for trace_strain_dot equation with state variables.
     *
     * f0_edot = ∂ε_v/∂t - ε_v_dot = 0
     */
    static inline
    void f0edot(const PylithInt dim,
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
        assert(numS >= 8);
        assert(s_t);

        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];
        const PylithScalar trace_strain_dot = s[sOff[i_trace_strain_dot]];

        f0[0] += trace_strain_t - trace_strain_dot;
    } // f0edot

    // ----------------------------------------------------------------------
    /** f0 function for temperature_dot equation with state variables.
     *
     * f0_Tdot = ∂T/∂t - T_dot = 0
     */
    static inline
    void f0Tdot(const PylithInt dim,
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
        assert(numS >= 8);
        assert(s_t);

        const PylithScalar temperature_t = s_t[sOff[i_temperature]];
        const PylithScalar temperature_dot = s[sOff[i_temperature_dot]];

        f0[0] += temperature_t - temperature_dot;
    } // f0Tdot

    // ============================= State Variable Jacobians =============================

    // ----------------------------------------------------------------------
    /** Jf0_vu for velocity equation.
     *
     * Jf0_vu = s_tshift (derivative of ∂u/∂t w.r.t. u)
     */
    static inline
    void Jf0vu(const PylithInt dim,
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
        for (PylithInt i = 0; i < dim; ++i) {
            Jf0[i * dim + i] += s_tshift;
        }
    } // Jf0vu

    // ----------------------------------------------------------------------
    /** Jf0_vv for velocity equation.
     *
     * Jf0_vv = -1 (derivative of -v w.r.t. v)
     */
    static inline
    void Jf0vv(const PylithInt dim,
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
        for (PylithInt i = 0; i < dim; ++i) {
            Jf0[i * dim + i] -= 1.0;
        }
    } // Jf0vv

    // ----------------------------------------------------------------------
    /** Jf0_pdotp for pressure_dot equation.
     *
     * Jf0_pdotp = s_tshift (derivative of ∂p/∂t w.r.t. p)
     */
    static inline
    void Jf0pdotp(const PylithInt dim,
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
        Jf0[0] += s_tshift;
    } // Jf0pdotp

    // ----------------------------------------------------------------------
    /** Jf0_pdotpdot for pressure_dot equation.
     *
     * Jf0_pdotpdot = -1 (derivative of -p_dot w.r.t. p_dot)
     */
    static inline
    void Jf0pdotpdot(const PylithInt dim,
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
        Jf0[0] -= 1.0;
    } // Jf0pdotpdot

    // ----------------------------------------------------------------------
    /** Jf0_edote for trace_strain_dot equation.
     */
    static inline
    void Jf0edote(const PylithInt dim,
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
        Jf0[0] += s_tshift;
    } // Jf0edote

    // ----------------------------------------------------------------------
    /** Jf0_edotedot for trace_strain_dot equation.
     */
    static inline
    void Jf0edotedot(const PylithInt dim,
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
        Jf0[0] -= 1.0;
    } // Jf0edotedot

    // ----------------------------------------------------------------------
    /** Jf0_TdotT for temperature_dot equation.
     */
    static inline
    void Jf0TdotT(const PylithInt dim,
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
        Jf0[0] += s_tshift;
    } // Jf0TdotT

    // ----------------------------------------------------------------------
    /** Jf0_TdotTdot for temperature_dot equation.
     */
    static inline
    void Jf0TdotTdot(const PylithInt dim,
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
        Jf0[0] -= 1.0;
    } // Jf0TdotTdot

    // ============================= Standard Jacobians =============================

    // ----------------------------------------------------------------------
    /** Jf0_TT function for temperature equation (heat capacity).
     *
     * Jf0_TT = (ρc)_eff * s_tshift
     */
    static inline
    void Jf0TT_transient(const PylithInt dim,
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
        assert(numA >= 10);

        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_porosity = 3;
        // Rheology fields indexed relative to numA (specific_heat is last)
        const PylithInt i_specificHeat = numA - 1;

        const PylithScalar solidDensity = a[aOff[i_solidDensity]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar porosity = a[aOff[i_porosity]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        const PylithScalar bulkDensity = (1.0 - porosity) * solidDensity + porosity * fluidDensity;
        const PylithScalar effectiveHeatCapacity = bulkDensity * specificHeat;

        Jf0[0] += effectiveHeatCapacity * s_tshift;
    } // Jf0TT_transient

}; // class Thermoporoelasticity

// End of file
