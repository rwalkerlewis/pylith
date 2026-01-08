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
 * Kernels for isotropic heat conduction.
 *
 * Heat flux (Fourier's law):
 * q = -k ∇T
 *
 * where:
 *   q = heat flux vector
 *   k = thermal conductivity (scalar for isotropic case)
 *   T = temperature
 *
 * In weak form, the heat flux term becomes:
 * ∫_Ω ∇v · k∇T dV
 *
 * f1 term: k * ∇T
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

class pylith::fekernels::IsotropicHeat {
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
    /** f1 function for isotropic heat conduction (LHS, implicit).
     *
     * f1 = k * ∇T
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, ...]
     */
    static inline
    void f1T(const PylithInt dim,
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
        // Solution field indices
        const PylithInt i_temperature = 0;

        // Auxiliary field indices
        const PylithInt i_thermalConductivity = 2;

        assert(numS >= 1);
        assert(numA >= 3);
        assert(sOff_x[i_temperature] >= 0);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        // f1 = k * ∇T
        for (PylithInt d = 0; d < dim; ++d) {
            f1[d] += thermalConductivity * temperature_x[d];
        } // for
    } // f1T

    // ============================= RHS Residual ================================

    // ----------------------------------------------------------------------
    /** g1 function for isotropic heat conduction (RHS, explicit).
     *
     * g1 = -k * ∇T
     *
     * Note: The negative sign is because in explicit formulation,
     * the heat flux term goes to the RHS with a sign change.
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, ...]
     */
    static inline
    void g1T(const PylithInt dim,
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
        // Solution field indices
        const PylithInt i_temperature = 0;

        // Auxiliary field indices
        const PylithInt i_thermalConductivity = 2;

        assert(numS >= 1);
        assert(numA >= 3);
        assert(sOff_x[i_temperature] >= 0);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        // g1 = -k * ∇T (note the negative sign for RHS)
        for (PylithInt d = 0; d < dim; ++d) {
            g1[d] -= thermalConductivity * temperature_x[d];
        } // for
    } // g1T

    // ============================= LHS Jacobian ================================

    // ----------------------------------------------------------------------
    /** Jf3 function for isotropic heat conduction (LHS Jacobian).
     *
     * Jf3_TT = k * δ_ij
     *
     * This is the Jacobian of f1 with respect to ∇T:
     * ∂(k * ∇T)/∂(∇T) = k * I (identity matrix)
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, ...]
     */
    static inline
    void Jf3TT(const PylithInt dim,
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
        // Auxiliary field indices
        const PylithInt i_thermalConductivity = 2;

        assert(numA >= 3);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        // Jf3 = k * δ_ij (diagonal entries only for isotropic case)
        // The Jacobian tensor is J[i][j] = ∂f1_i/∂(∇T)_j = k * δ_ij
        for (PylithInt d = 0; d < dim; ++d) {
            Jf3[d * dim + d] += thermalConductivity;
        } // for
    } // Jf3TT

    // ============================= Derived Fields ================================

    // ----------------------------------------------------------------------
    /** Compute heat flux as a vector field for output.
     *
     * q = -k * ∇T
     *
     * Solution fields: [temperature]
     * Auxiliary fields: [density, specific_heat, thermal_conductivity, ...]
     */
    static inline
    void heatFlux_asVector(const PylithInt dim,
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
                           PylithScalar heatFlux[]) {
        // Solution field indices
        const PylithInt i_temperature = 0;

        // Auxiliary field indices
        const PylithInt i_thermalConductivity = 2;

        assert(numS >= 1);
        assert(numA >= 3);
        assert(sOff_x[i_temperature] >= 0);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        // q = -k * ∇T (Fourier's law)
        for (PylithInt d = 0; d < dim; ++d) {
            heatFlux[d] = -thermalConductivity * temperature_x[d];
        } // for
    } // heatFlux_asVector

}; // class IsotropicHeat

// End of file
