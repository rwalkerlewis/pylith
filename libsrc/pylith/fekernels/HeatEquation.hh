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
 * Kernels for the heat equation.
 *
 * Heat Equation:
 * ρc ∂T/∂t - ∇·(k∇T) = Q
 *
 * where:
 *   T = temperature
 *   ρ = density
 *   c = specific heat capacity
 *   k = thermal conductivity
 *   Q = heat source
 *
 * Weak form (quasistatic):
 * ∫_Ω ∇v · k∇T dV = ∫_Ω v Q dV
 *
 * Weak form (time-dependent):
 * ∫_Ω v ρc ∂T/∂t dV + ∫_Ω ∇v · k∇T dV = ∫_Ω v Q dV
 *
 * In PyLith's formulation (F(t,s,\dot{s}) = G(t)):
 * LHS: f0 = ρc T_t, f1 = k∇T
 * RHS: g0 = Q, g1 = 0
 *
 *** Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0/f1 [numComponents] Output residual.
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::HeatEquation {
public:

    // ============================= Solution Fields =============================
    // Solution fields: [temperature]
    //
    // Auxiliary fields:
    //   Index 0: density
    //   Index 1: specific_heat
    //   Index 2: thermal_conductivity
    //   Index 3 (optional): heat_source

    // ============================= LHS Residual ================================

    // ----------------------------------------------------------------------
    /** f0 function for heat equation, time-dependent (LHS).
     *
     * f0 = ρc T_t
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, ...]
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
        // Solution field indices
        const PylithInt i_temperature = 0;

        // Auxiliary field indices
        const PylithInt i_density = 0;
        const PylithInt i_specificHeat = 1;

        assert(numS >= 1);
        assert(numA >= 3);
        assert(sOff[i_temperature] >= 0);
        assert(aOff[i_density] >= 0);
        assert(aOff[i_specificHeat] >= 0);
        assert(s_t);

        const PylithScalar temperature_t = s_t[sOff[i_temperature]];
        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        f0[0] += density * specificHeat * temperature_t;
    } // f0T_timedep

    // ----------------------------------------------------------------------
    /** f0 function for heat equation with heat source (LHS, quasistatic).
     *
     * f0 = -Q (moved to LHS)
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, heat_source]
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
        // Auxiliary field indices
        const PylithInt i_heatSource = 3;

        assert(numA >= 4);
        assert(aOff[i_heatSource] >= 0);

        const PylithScalar heatSource = a[aOff[i_heatSource]];

        // Negative sign because we move to LHS: F(t,s) = 0 => -Q on LHS
        f0[0] -= heatSource;
    } // f0T_source

    // ============================= RHS Residual ================================

    // ----------------------------------------------------------------------
    /** g0 function for heat equation with heat source (RHS).
     *
     * g0 = Q
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, heat_source]
     */
    static inline
    void g0T_source(const PylithInt dim,
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
        // Auxiliary field indices
        const PylithInt i_heatSource = 3;

        assert(numA >= 4);
        assert(aOff[i_heatSource] >= 0);

        const PylithScalar heatSource = a[aOff[i_heatSource]];

        g0[0] += heatSource;
    } // g0T_source

    // ============================= LHS Jacobian ================================

    // ----------------------------------------------------------------------
    /** Jf0 function for heat equation, time-dependent (LHS).
     *
     * Jf0_TT = ρc * s_tshift
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity]
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
        // Auxiliary field indices
        const PylithInt i_density = 0;
        const PylithInt i_specificHeat = 1;

        assert(numA >= 3);
        assert(aOff[i_density] >= 0);
        assert(aOff[i_specificHeat] >= 0);

        const PylithScalar density = a[aOff[i_density]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        Jf0[0] += density * specificHeat * s_tshift;
    } // Jf0TT_timedep

}; // class HeatEquation

// End of file
