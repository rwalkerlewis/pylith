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
 * Kernels for isotropic linear thermoelasticity.
 *
 * Constitutive relation with thermal strain (plane strain / 3D):
 * σ = λ tr(ε - ε_th) I + 2μ (ε - ε_th)
 * where ε_th = α (T - T_ref) I
 *
 * So:
 * σ = λ (tr(ε) - 3α(T - T_ref)) I + 2μ (ε - α(T - T_ref)I)
 *   = λ tr(ε) I + 2μ ε - (3λ + 2μ) α (T - T_ref) I
 *   = λ tr(ε) I + 2μ ε - 3K α (T - T_ref) I
 * where K = λ + 2μ/3 is the bulk modulus
 *
 * Heat flux (Fourier's law):
 * q = -k ∇T
 *
 * Solution fields: [displacement, temperature]
 * Auxiliary fields (base):
 *   0: density
 *   1: specific_heat
 *   2: thermal_conductivity
 *   3: reference_temperature
 *   4: thermal_expansion_coefficient
 *   5: shear_modulus
 *   6: bulk_modulus
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert>

class pylith::fekernels::IsotropicLinearThermoelasticity {
public:

    // ============================= Auxiliary Field Indices =============================
    // Base fields
    static const PylithInt i_density = 0;
    static const PylithInt i_specificHeat = 1;
    static const PylithInt i_thermalConductivity = 2;
    static const PylithInt i_refTemperature = 3;
    static const PylithInt i_thermalExpansionCoeff = 4;
    static const PylithInt i_shearModulus = 5;
    static const PylithInt i_bulkModulus = 6;

    // Solution field indices
    static const PylithInt i_displacement = 0;
    static const PylithInt i_temperature = 1;

    // ============================= Displacement Equation: Stress ================================

    // ----------------------------------------------------------------------
    /** f1 function for stress in displacement equation (LHS, implicit).
     *
     * f1_u = σ = λ tr(ε) I + 2μ ε - 3K α (T - T_ref) I
     *
     * Plane strain formulation.
     */
    static inline
    void f1u_PlaneStrain(const PylithInt dim,
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
        assert(2 == dim);
        assert(numS >= 2);
        assert(numA >= 7);
        assert(sOff_x[i_displacement] >= 0);
        assert(sOff[i_temperature] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_refTemperature] >= 0);
        assert(aOff[i_thermalExpansionCoeff] >= 0);

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar refTemperature = a[aOff[i_refTemperature]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpansionCoeff]];

        // Strain components
        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy;

        // Thermal strain contribution
        const PylithScalar deltaT = temperature - refTemperature;
        const PylithScalar thermalStress = 3.0 * bulkModulus * thermalExpansionCoeff * deltaT;

        // Lame parameters: λ = K - 2μ/3
        const PylithScalar lambda = bulkModulus - 2.0 * shearModulus / 3.0;

        // Stress (plane strain): σ_zz ≠ 0 but ε_zz = 0
        // σ = λ tr(ε) I + 2μ ε - thermalStress * I
        const PylithScalar stress_xx = lambda * trace_strain + 2.0 * shearModulus * strain_xx - thermalStress;
        const PylithScalar stress_yy = lambda * trace_strain + 2.0 * shearModulus * strain_yy - thermalStress;
        const PylithScalar stress_xy = 2.0 * shearModulus * strain_xy;

        // f1 is stored as a 2x2 tensor (row-major)
        f1[0*dim+0] += stress_xx;
        f1[0*dim+1] += stress_xy;
        f1[1*dim+0] += stress_xy;
        f1[1*dim+1] += stress_yy;
    } // f1u_PlaneStrain

    // ----------------------------------------------------------------------
    /** f1 function for stress in displacement equation (LHS, implicit) - 3D.
     *
     * f1_u = σ = λ tr(ε) I + 2μ ε - 3K α (T - T_ref) I
     */
    static inline
    void f1u_3D(const PylithInt dim,
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
        assert(3 == dim);
        assert(numS >= 2);
        assert(numA >= 7);
        assert(sOff_x[i_displacement] >= 0);
        assert(sOff[i_temperature] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_refTemperature] >= 0);
        assert(aOff[i_thermalExpansionCoeff] >= 0);

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar refTemperature = a[aOff[i_refTemperature]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpansionCoeff]];

        // Strain components
        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_zz = disp_x[2*dim+2];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar strain_yz = 0.5 * (disp_x[1*dim+2] + disp_x[2*dim+1]);
        const PylithScalar strain_xz = 0.5 * (disp_x[0*dim+2] + disp_x[2*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy + strain_zz;

        // Thermal strain contribution
        const PylithScalar deltaT = temperature - refTemperature;
        const PylithScalar thermalStress = 3.0 * bulkModulus * thermalExpansionCoeff * deltaT;

        // Lame parameters
        const PylithScalar lambda = bulkModulus - 2.0 * shearModulus / 3.0;

        // Stress
        const PylithScalar stress_xx = lambda * trace_strain + 2.0 * shearModulus * strain_xx - thermalStress;
        const PylithScalar stress_yy = lambda * trace_strain + 2.0 * shearModulus * strain_yy - thermalStress;
        const PylithScalar stress_zz = lambda * trace_strain + 2.0 * shearModulus * strain_zz - thermalStress;
        const PylithScalar stress_xy = 2.0 * shearModulus * strain_xy;
        const PylithScalar stress_yz = 2.0 * shearModulus * strain_yz;
        const PylithScalar stress_xz = 2.0 * shearModulus * strain_xz;

        // f1 stored as 3x3 tensor (row-major)
        f1[0*dim+0] += stress_xx;
        f1[0*dim+1] += stress_xy;
        f1[0*dim+2] += stress_xz;
        f1[1*dim+0] += stress_xy;
        f1[1*dim+1] += stress_yy;
        f1[1*dim+2] += stress_yz;
        f1[2*dim+0] += stress_xz;
        f1[2*dim+1] += stress_yz;
        f1[2*dim+2] += stress_zz;
    } // f1u_3D

    // ============================= Temperature Equation: Heat Flux ================================

    // ----------------------------------------------------------------------
    /** f1 function for heat flux in temperature equation (LHS, implicit).
     *
     * f1_T = k ∇T
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
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff_x[i_temperature] >= 0);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        for (PylithInt d = 0; d < dim; ++d) {
            f1[d] += thermalConductivity * temperature_x[d];
        } // for
    } // f1T

    // ============================= LHS Jacobian ================================

    // ----------------------------------------------------------------------
    /** Jf3_uu function for elastic stiffness (plane strain).
     *
     * Jf3_uu = C_ijkl (elastic stiffness tensor)
     */
    static inline
    void Jf3uu_PlaneStrain(const PylithInt dim,
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
        assert(2 == dim);
        assert(numA >= 7);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar lambda = bulkModulus - 2.0 * shearModulus / 3.0;

        // Stiffness matrix for plane strain (4x4 in Voigt notation mapped to tensor)
        // C_ijkl: Jf3[((i*dim+j)*dim+k)*dim+l]
        const PylithReal C1111 = lambda + 2.0 * shearModulus;
        const PylithReal C1122 = lambda;
        const PylithReal C1212 = shearModulus;
        const PylithReal C2211 = lambda;
        const PylithReal C2222 = lambda + 2.0 * shearModulus;
        const PylithReal C2121 = shearModulus;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         * 0: j0000 = C1111
         * 1: j0001 = C1112 = 0
         * 2: j0010 = C1211 = 0
         * 3: j0011 = C1212
         * 4: j0100 = C1121 = 0
         * 5: j0101 = C1122
         * 6: j0110 = C1221 = shearModulus
         * 7: j0111 = C1222 = 0
         * 8: j1000 = C2111 = 0
         * 9: j1001 = C2112 = shearModulus
         * 10: j1010 = C2211
         * 11: j1011 = C2212 = 0
         * 12: j1100 = C2121
         * 13: j1101 = C2122 = 0
         * 14: j1110 = C2221 = 0
         * 15: j1111 = C2222
         */
        Jf3[0] += C1111; // j0000
        Jf3[3] += C1212; // j0011
        Jf3[5] += C1122; // j0101
        Jf3[6] += C1212; // j0110
        Jf3[9] += C2121; // j1001
        Jf3[10] += C2211; // j1010
        Jf3[12] += C2121; // j1100
        Jf3[15] += C2222; // j1111
    } // Jf3uu_PlaneStrain

    // ----------------------------------------------------------------------
    /** Jf3_uu function for elastic stiffness (3D).
     */
    static inline
    void Jf3uu_3D(const PylithInt dim,
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
        assert(3 == dim);
        assert(numA >= 7);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_bulkModulus] >= 0);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar lambda = bulkModulus - 2.0 * shearModulus / 3.0;

        const PylithReal C1111 = lambda + 2.0 * shearModulus;
        const PylithReal C1122 = lambda;
        const PylithReal C1133 = lambda;
        const PylithReal C1212 = shearModulus;
        const PylithReal C1313 = shearModulus;
        const PylithReal C2222 = lambda + 2.0 * shearModulus;
        const PylithReal C2233 = lambda;
        const PylithReal C2323 = shearModulus;
        const PylithReal C3333 = lambda + 2.0 * shearModulus;

        // Full 3D stiffness tensor (81 components, but only non-zero ones matter)
        // Jf3[((i*dim+j)*dim+k)*dim+l] = C_ijkl
        Jf3[0] += C1111;   // 0000
        Jf3[4] += C1212;   // 0011
        Jf3[8] += C1313;   // 0022
        Jf3[10] += C1122;  // 0101
        Jf3[12] += C1212;  // 0110
        Jf3[20] += C1133;  // 0202
        Jf3[24] += C1313;  // 0220
        Jf3[28] += C1212;  // 1001
        Jf3[30] += C1122;  // 1010
        Jf3[40] += C2222;  // 1111
        Jf3[44] += C2323;  // 1122
        Jf3[50] += C2233;  // 1212
        Jf3[52] += C2323;  // 1221
        Jf3[56] += C1313;  // 2002
        Jf3[60] += C1133;  // 2020
        Jf3[68] += C2323;  // 2112
        Jf3[70] += C2233;  // 2121
        Jf3[80] += C3333;  // 2222
    } // Jf3uu_3D

    // ----------------------------------------------------------------------
    /** Jf2_uT function for coupling: stress depends on temperature.
     *
     * Jf2_uT = -3K α I (derivative of stress w.r.t. temperature)
     *
     * Plane strain version.
     */
    static inline
    void Jf2uT_PlaneStrain(const PylithInt dim,
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
                           PylithScalar Jf2[]) {
        assert(2 == dim);
        assert(numA >= 7);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_thermalExpansionCoeff] >= 0);

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpansionCoeff]];

        // ∂σ/∂T = -3K α I (negative because σ = ... - 3Kα(T-Tref)I)
        // Jf2[i*dim+j] = ∂f1_ij/∂T
        const PylithScalar coupling = -3.0 * bulkModulus * thermalExpansionCoeff;

        Jf2[0*dim+0] += coupling; // ∂σ_xx/∂T
        Jf2[1*dim+1] += coupling; // ∂σ_yy/∂T
    } // Jf2uT_PlaneStrain

    // ----------------------------------------------------------------------
    /** Jf2_uT function for 3D.
     */
    static inline
    void Jf2uT_3D(const PylithInt dim,
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
                  PylithScalar Jf2[]) {
        assert(3 == dim);
        assert(numA >= 7);
        assert(aOff[i_bulkModulus] >= 0);
        assert(aOff[i_thermalExpansionCoeff] >= 0);

        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpansionCoeff]];

        const PylithScalar coupling = -3.0 * bulkModulus * thermalExpansionCoeff;

        Jf2[0*dim+0] += coupling; // ∂σ_xx/∂T
        Jf2[1*dim+1] += coupling; // ∂σ_yy/∂T
        Jf2[2*dim+2] += coupling; // ∂σ_zz/∂T
    } // Jf2uT_3D

    // ----------------------------------------------------------------------
    /** Jf3_TT function for thermal conductivity.
     *
     * Jf3_TT = k I
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
        assert(numA >= 3);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        for (PylithInt d = 0; d < dim; ++d) {
            Jf3[d * dim + d] += thermalConductivity;
        } // for
    } // Jf3TT

    // ============================= Derived Fields ================================

    // ----------------------------------------------------------------------
    /** Compute Cauchy stress as a vector for output (plane strain).
     */
    static inline
    void cauchyStress_PlaneStrain(const PylithInt dim,
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
                                  PylithScalar stress[]) {
        assert(2 == dim);
        assert(numS >= 2);
        assert(numA >= 7);

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
        const PylithScalar refTemperature = a[aOff[i_refTemperature]];
        const PylithScalar thermalExpansionCoeff = a[aOff[i_thermalExpansionCoeff]];

        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy;

        const PylithScalar deltaT = temperature - refTemperature;
        const PylithScalar thermalStress = 3.0 * bulkModulus * thermalExpansionCoeff * deltaT;
        const PylithScalar lambda = bulkModulus - 2.0 * shearModulus / 3.0;

        // Output in Voigt notation: [σ_xx, σ_yy, σ_zz, σ_xy]
        stress[0] = lambda * trace_strain + 2.0 * shearModulus * strain_xx - thermalStress; // σ_xx
        stress[1] = lambda * trace_strain + 2.0 * shearModulus * strain_yy - thermalStress; // σ_yy
        stress[2] = lambda * trace_strain - thermalStress; // σ_zz (plane strain: ε_zz = 0)
        stress[3] = 2.0 * shearModulus * strain_xy; // σ_xy
    } // cauchyStress_PlaneStrain

    // ----------------------------------------------------------------------
    /** Compute heat flux as a vector for output.
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
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff_x[i_temperature] >= 0);
        assert(aOff[i_thermalConductivity] >= 0);

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalConductivity]];

        // q = -k ∇T
        for (PylithInt d = 0; d < dim; ++d) {
            heatFlux[d] = -thermalConductivity * temperature_x[d];
        } // for
    } // heatFlux_asVector

}; // class IsotropicLinearThermoelasticity

// End of file
