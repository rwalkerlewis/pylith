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
 * Kernels for isotropic linear thermoporoelasticity.
 *
 * Constitutive equations:
 *
 * 1. Stress-strain relation with thermal and pore pressure effects:
 *    σ = λ_d tr(ε) I + 2μ ε - 3K_d α_T (T - T_ref) I - α_B p I
 *    where:
 *    - λ_d, μ = drained Lame parameters
 *    - K_d = drained bulk modulus
 *    - α_T = linear thermal expansion coefficient
 *    - α_B = Biot coefficient
 *    - p = pore pressure
 *
 * 2. Fluid content:
 *    ζ = α_B ε_v + p/M + 3α_f (T - T_ref)
 *    where:
 *    - ε_v = volumetric strain = tr(ε)
 *    - M = Biot modulus
 *    - α_f = fluid thermal expansivity × porosity
 *
 * 3. Darcy flux:
 *    q = -k/μ_f (∇p - ρ_f g)
 *
 * 4. Heat flux:
 *    q_T = -k_T ∇T
 *
 * Solution fields: [displacement, pressure, trace_strain, temperature]
 *
 * Auxiliary fields (indices may vary based on optional fields):
 *   0: solid_density
 *   1: fluid_density
 *   2: fluid_viscosity
 *   3: porosity
 *   (optional): gravity_field, body_force, source_density, heat_source
 *   Rheology: biot_coefficient, biot_modulus, drained_bulk_modulus, shear_modulus,
 *             reference_temperature, thermal_expansion_coefficient, fluid_thermal_expansion,
 *             thermal_conductivity, specific_heat, permeability
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Tensor.hh"

#include "pylith/utils/types.hh"

#include <cassert>

class pylith::fekernels::IsotropicLinearThermoporoelasticity {
public:

    // ============================= Auxiliary Field Indices =============================
    // Base poroelastic fields
    static const PylithInt i_solidDensity = 0;
    static const PylithInt i_fluidDensity = 1;
    static const PylithInt i_fluidViscosity = 2;
    static const PylithInt i_porosity = 3;

    // Solution field indices
    static const PylithInt i_displacement = 0;
    static const PylithInt i_pressure = 1;
    static const PylithInt i_trace_strain = 2;
    static const PylithInt i_temperature = 3;

    // ============================= Stress Kernels =============================

    // ----------------------------------------------------------------------
    /** f1 function for stress in displacement equation (LHS, implicit) - Plane Strain.
     *
     * σ = λ_d tr(ε) I + 2μ ε - 3K_d α_T (T - T_ref) I - α_B p I
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
        assert(numS >= 4);

        // Constants array provides auxiliary field indices for rheology-specific fields
        // constants[0]: i_biotCoeff, constants[1]: i_drainedBulkMod, constants[2]: i_shearMod
        // constants[3]: i_refTemp, constants[4]: i_thermalExpCoeff
        const PylithInt i_biotCoeff = numConstants >= 5 ? (PylithInt)constants[0] : 4;
        const PylithInt i_drainedBulkMod = numConstants >= 5 ? (PylithInt)constants[1] : 5;
        const PylithInt i_shearMod = numConstants >= 5 ? (PylithInt)constants[2] : 6;
        const PylithInt i_refTemp = numConstants >= 5 ? (PylithInt)constants[3] : 7;
        const PylithInt i_thermalExpCoeff = numConstants >= 5 ? (PylithInt)constants[4] : 8;

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];
        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar shearMod = a[aOff[i_shearMod]];
        const PylithScalar refTemp = a[aOff[i_refTemp]];
        const PylithScalar thermalExpCoeff = a[aOff[i_thermalExpCoeff]];

        // Strain components
        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy;

        // Drained Lame parameter
        const PylithScalar lambda_d = drainedBulkMod - 2.0 * shearMod / 3.0;

        // Temperature effect
        const PylithScalar deltaT = temperature - refTemp;
        const PylithScalar thermalStress = 3.0 * drainedBulkMod * thermalExpCoeff * deltaT;

        // Stress components: σ = λ_d tr(ε) I + 2μ ε - thermalStress I - α_B p I
        const PylithScalar stress_xx = lambda_d * trace_strain + 2.0 * shearMod * strain_xx 
                                     - thermalStress - biotCoeff * pressure;
        const PylithScalar stress_yy = lambda_d * trace_strain + 2.0 * shearMod * strain_yy 
                                     - thermalStress - biotCoeff * pressure;
        const PylithScalar stress_xy = 2.0 * shearMod * strain_xy;

        f1[0*dim+0] += stress_xx;
        f1[0*dim+1] += stress_xy;
        f1[1*dim+0] += stress_xy;
        f1[1*dim+1] += stress_yy;
    } // f1u_PlaneStrain

    // ----------------------------------------------------------------------
    /** f1 function for stress in displacement equation (LHS, implicit) - 3D.
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
        assert(numS >= 4);

        const PylithInt i_biotCoeff = numConstants >= 5 ? (PylithInt)constants[0] : 4;
        const PylithInt i_drainedBulkMod = numConstants >= 5 ? (PylithInt)constants[1] : 5;
        const PylithInt i_shearMod = numConstants >= 5 ? (PylithInt)constants[2] : 6;
        const PylithInt i_refTemp = numConstants >= 5 ? (PylithInt)constants[3] : 7;
        const PylithInt i_thermalExpCoeff = numConstants >= 5 ? (PylithInt)constants[4] : 8;

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];
        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar shearMod = a[aOff[i_shearMod]];
        const PylithScalar refTemp = a[aOff[i_refTemp]];
        const PylithScalar thermalExpCoeff = a[aOff[i_thermalExpCoeff]];

        // Strain components
        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_zz = disp_x[2*dim+2];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar strain_yz = 0.5 * (disp_x[1*dim+2] + disp_x[2*dim+1]);
        const PylithScalar strain_xz = 0.5 * (disp_x[0*dim+2] + disp_x[2*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy + strain_zz;

        const PylithScalar lambda_d = drainedBulkMod - 2.0 * shearMod / 3.0;
        const PylithScalar deltaT = temperature - refTemp;
        const PylithScalar thermalStress = 3.0 * drainedBulkMod * thermalExpCoeff * deltaT;

        const PylithScalar stress_xx = lambda_d * trace_strain + 2.0 * shearMod * strain_xx 
                                     - thermalStress - biotCoeff * pressure;
        const PylithScalar stress_yy = lambda_d * trace_strain + 2.0 * shearMod * strain_yy 
                                     - thermalStress - biotCoeff * pressure;
        const PylithScalar stress_zz = lambda_d * trace_strain + 2.0 * shearMod * strain_zz 
                                     - thermalStress - biotCoeff * pressure;
        const PylithScalar stress_xy = 2.0 * shearMod * strain_xy;
        const PylithScalar stress_yz = 2.0 * shearMod * strain_yz;
        const PylithScalar stress_xz = 2.0 * shearMod * strain_xz;

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

    // ============================= Fluid Content / Pressure Kernels =============================

    // ----------------------------------------------------------------------
    /** f0 function for pressure equation (fluid content time derivative).
     *
     * f0_p = ∂ζ/∂t = α_B ∂ε_v/∂t + (1/M) ∂p/∂t + 3α_f ∂T/∂t
     *
     * Using trace_strain as ε_v and its time derivative.
     */
    static inline
    void f0p_fluidContent(const PylithInt dim,
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
        assert(numS >= 4);
        assert(s_t);

        // Constants array provides indices
        const PylithInt i_biotCoeff = numConstants >= 3 ? (PylithInt)constants[0] : 4;
        const PylithInt i_biotMod = numConstants >= 3 ? (PylithInt)constants[1] : 5;
        const PylithInt i_fluidThermalExp = numConstants >= 3 ? (PylithInt)constants[2] : 9;

        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];
        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar temperature_t = s_t[sOff[i_temperature]];

        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];
        const PylithScalar biotMod = a[aOff[i_biotMod]];
        const PylithScalar fluidThermalExp = a[aOff[i_fluidThermalExp]];

        // ∂ζ/∂t = α_B ∂ε_v/∂t + (1/M) ∂p/∂t + 3α_f ∂T/∂t
        f0[0] += biotCoeff * trace_strain_t + pressure_t / biotMod + 3.0 * fluidThermalExp * temperature_t;
    } // f0p_fluidContent

    // ----------------------------------------------------------------------
    /** f0 function for pressure equation with source density.
     */
    static inline
    void f0p_source(const PylithInt dim,
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
        // Add fluid content term
        f0p_fluidContent(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0);

        // Subtract source density (it's on RHS)
        const PylithInt i_sourceDensity = numConstants >= 4 ? (PylithInt)constants[3] : 10;
        const PylithScalar sourceDensity = a[aOff[i_sourceDensity]];

        f0[0] -= sourceDensity;
    } // f0p_source

    // ----------------------------------------------------------------------
    /** f1 function for Darcy flux in pressure equation.
     *
     * f1_p = q = -k/μ (∇p - ρ_f g)
     *
     * Without gravity.
     */
    static inline
    void f1p_darcy(const PylithInt dim,
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
        const PylithInt i_permeability = numConstants >= 2 ? (PylithInt)constants[0] : 10;
        const PylithInt i_fluidVisc = i_fluidViscosity;

        const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar permeability = a[aOff[i_permeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidVisc]];

        // Darcy flux: q = -k/μ ∇p
        const PylithScalar darcyCoeff = permeability / fluidViscosity;

        for (PylithInt d = 0; d < dim; ++d) {
            f1[d] -= darcyCoeff * pressure_x[d];
        }
    } // f1p_darcy

    // ----------------------------------------------------------------------
    /** f1 function for Darcy flux with gravity.
     */
    static inline
    void f1p_darcy_grav(const PylithInt dim,
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
        const PylithInt i_permeability = numConstants >= 3 ? (PylithInt)constants[0] : 10;
        const PylithInt i_gravity = numConstants >= 3 ? (PylithInt)constants[1] : 11;
        const PylithInt i_fluidVisc = i_fluidViscosity;
        const PylithInt i_fluidDens = i_fluidDensity;

        const PylithScalar* pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar permeability = a[aOff[i_permeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidVisc]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDens]];
        const PylithScalar* gravityField = &a[aOff[i_gravity]];

        const PylithScalar darcyCoeff = permeability / fluidViscosity;

        // q = -k/μ (∇p - ρ_f g)
        for (PylithInt d = 0; d < dim; ++d) {
            f1[d] -= darcyCoeff * (pressure_x[d] - fluidDensity * gravityField[d]);
        }
    } // f1p_darcy_grav

    // ============================= Temperature Equation =============================

    // ----------------------------------------------------------------------
    /** f1 function for heat flux in temperature equation.
     *
     * f1_T = k_T ∇T
     */
    static inline
    void f1T_heatflux(const PylithInt dim,
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
        const PylithInt i_thermalCond = numConstants >= 1 ? (PylithInt)constants[0] : 11;

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalCond]];

        for (PylithInt d = 0; d < dim; ++d) {
            f1[d] += thermalConductivity * temperature_x[d];
        }
    } // f1T_heatflux

    // ============================= Jacobians =============================

    // ----------------------------------------------------------------------
    /** Jf3_uu for elastic stiffness (plane strain).
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

        const PylithInt i_drainedBulkMod = numConstants >= 2 ? (PylithInt)constants[0] : 5;
        const PylithInt i_shearMod = numConstants >= 2 ? (PylithInt)constants[1] : 6;

        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar shearMod = a[aOff[i_shearMod]];
        const PylithScalar lambda_d = drainedBulkMod - 2.0 * shearMod / 3.0;

        const PylithReal C1111 = lambda_d + 2.0 * shearMod;
        const PylithReal C1122 = lambda_d;
        const PylithReal C1212 = shearMod;
        const PylithReal C2222 = lambda_d + 2.0 * shearMod;

        Jf3[0] += C1111;  // j0000
        Jf3[3] += C1212;  // j0011
        Jf3[5] += C1122;  // j0101
        Jf3[6] += C1212;  // j0110
        Jf3[9] += C1212;  // j1001
        Jf3[10] += C1122; // j1010
        Jf3[12] += C1212; // j1100
        Jf3[15] += C2222; // j1111
    } // Jf3uu_PlaneStrain

    // ----------------------------------------------------------------------
    /** Jf3_uu for elastic stiffness (3D).
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

        const PylithInt i_drainedBulkMod = numConstants >= 2 ? (PylithInt)constants[0] : 5;
        const PylithInt i_shearMod = numConstants >= 2 ? (PylithInt)constants[1] : 6;

        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar shearMod = a[aOff[i_shearMod]];
        const PylithScalar lambda_d = drainedBulkMod - 2.0 * shearMod / 3.0;

        const PylithReal C1111 = lambda_d + 2.0 * shearMod;
        const PylithReal C1122 = lambda_d;
        const PylithReal C1212 = shearMod;

        Jf3[0] += C1111;   // 0000
        Jf3[4] += C1212;   // 0011
        Jf3[8] += C1212;   // 0022
        Jf3[10] += C1122;  // 0101
        Jf3[12] += C1212;  // 0110
        Jf3[20] += C1122;  // 0202
        Jf3[24] += C1212;  // 0220
        Jf3[28] += C1212;  // 1001
        Jf3[30] += C1122;  // 1010
        Jf3[40] += C1111;  // 1111
        Jf3[44] += C1212;  // 1122
        Jf3[50] += C1122;  // 1212
        Jf3[52] += C1212;  // 1221
        Jf3[56] += C1212;  // 2002
        Jf3[60] += C1122;  // 2020
        Jf3[68] += C1212;  // 2112
        Jf3[70] += C1122;  // 2121
        Jf3[80] += C1111;  // 2222
    } // Jf3uu_3D

    // ----------------------------------------------------------------------
    /** Jf2_up for Biot coupling (stress depends on pressure).
     *
     * ∂σ/∂p = -α_B I
     */
    static inline
    void Jf2up(const PylithInt dim,
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
        const PylithInt i_biotCoeff = numConstants >= 1 ? (PylithInt)constants[0] : 4;
        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];

        for (PylithInt d = 0; d < dim; ++d) {
            Jf2[d * dim + d] -= biotCoeff;
        }
    } // Jf2up

    // ----------------------------------------------------------------------
    /** Jf2_uT for thermal coupling (stress depends on temperature).
     *
     * ∂σ/∂T = -3K_d α_T I
     */
    static inline
    void Jf2uT(const PylithInt dim,
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
        const PylithInt i_drainedBulkMod = numConstants >= 2 ? (PylithInt)constants[0] : 5;
        const PylithInt i_thermalExpCoeff = numConstants >= 2 ? (PylithInt)constants[1] : 8;

        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar thermalExpCoeff = a[aOff[i_thermalExpCoeff]];

        const PylithScalar coupling = -3.0 * drainedBulkMod * thermalExpCoeff;

        for (PylithInt d = 0; d < dim; ++d) {
            Jf2[d * dim + d] += coupling;
        }
    } // Jf2uT

    // ----------------------------------------------------------------------
    /** Jf0_pp for pressure equation (storage coefficient).
     *
     * ∂f0_p/∂p_t = 1/M
     */
    static inline
    void Jf0pp(const PylithInt dim,
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
        const PylithInt i_biotMod = numConstants >= 1 ? (PylithInt)constants[0] : 5;
        const PylithScalar biotMod = a[aOff[i_biotMod]];

        Jf0[0] += s_tshift / biotMod;
    } // Jf0pp

    // ----------------------------------------------------------------------
    /** Jf0_pe for pressure-trace_strain coupling.
     *
     * ∂f0_p/∂ε_v_t = α_B
     */
    static inline
    void Jf0pe(const PylithInt dim,
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
        const PylithInt i_biotCoeff = numConstants >= 1 ? (PylithInt)constants[0] : 4;
        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];

        Jf0[0] += biotCoeff * s_tshift;
    } // Jf0pe

    // ----------------------------------------------------------------------
    /** Jf0_pT for pressure-temperature coupling.
     *
     * ∂f0_p/∂T_t = 3α_f
     */
    static inline
    void Jf0pT(const PylithInt dim,
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
        const PylithInt i_fluidThermalExp = numConstants >= 1 ? (PylithInt)constants[0] : 9;
        const PylithScalar fluidThermalExp = a[aOff[i_fluidThermalExp]];

        Jf0[0] += 3.0 * fluidThermalExp * s_tshift;
    } // Jf0pT

    // ----------------------------------------------------------------------
    /** Jf3_pp for Darcy conductivity.
     *
     * Jf3_pp = k/μ I
     */
    static inline
    void Jf3pp(const PylithInt dim,
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
        const PylithInt i_permeability = numConstants >= 1 ? (PylithInt)constants[0] : 10;
        const PylithScalar permeability = a[aOff[i_permeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar darcyCoeff = permeability / fluidViscosity;

        for (PylithInt d = 0; d < dim; ++d) {
            Jf3[d * dim + d] -= darcyCoeff;
        }
    } // Jf3pp

    // ----------------------------------------------------------------------
    /** Jf0_TT for temperature equation (heat capacity).
     */
    static inline
    void Jf0TT(const PylithInt dim,
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
        const PylithInt i_specificHeat = numConstants >= 1 ? (PylithInt)constants[0] : 12;

        const PylithScalar solidDensity = a[aOff[i_solidDensity]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar porosity = a[aOff[i_porosity]];
        const PylithScalar specificHeat = a[aOff[i_specificHeat]];

        const PylithScalar bulkDensity = (1.0 - porosity) * solidDensity + porosity * fluidDensity;
        const PylithScalar effectiveHeatCapacity = bulkDensity * specificHeat;

        Jf0[0] += effectiveHeatCapacity * s_tshift;
    } // Jf0TT

    // ----------------------------------------------------------------------
    /** Jf3_TT for thermal conductivity.
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
        const PylithInt i_thermalCond = numConstants >= 1 ? (PylithInt)constants[0] : 11;
        const PylithScalar thermalConductivity = a[aOff[i_thermalCond]];

        for (PylithInt d = 0; d < dim; ++d) {
            Jf3[d * dim + d] += thermalConductivity;
        }
    } // Jf3TT

    // ============================= Derived Fields =============================

    // ----------------------------------------------------------------------
    /** Compute Cauchy stress (plane strain).
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

        const PylithInt i_biotCoeff = numConstants >= 5 ? (PylithInt)constants[0] : 4;
        const PylithInt i_drainedBulkMod = numConstants >= 5 ? (PylithInt)constants[1] : 5;
        const PylithInt i_shearMod = numConstants >= 5 ? (PylithInt)constants[2] : 6;
        const PylithInt i_refTemp = numConstants >= 5 ? (PylithInt)constants[3] : 7;
        const PylithInt i_thermalExpCoeff = numConstants >= 5 ? (PylithInt)constants[4] : 8;

        const PylithScalar* disp_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];
        const PylithScalar drainedBulkMod = a[aOff[i_drainedBulkMod]];
        const PylithScalar shearMod = a[aOff[i_shearMod]];
        const PylithScalar refTemp = a[aOff[i_refTemp]];
        const PylithScalar thermalExpCoeff = a[aOff[i_thermalExpCoeff]];

        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_xy = 0.5 * (disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar trace_strain = strain_xx + strain_yy;

        const PylithScalar lambda_d = drainedBulkMod - 2.0 * shearMod / 3.0;
        const PylithScalar deltaT = temperature - refTemp;
        const PylithScalar thermalStress = 3.0 * drainedBulkMod * thermalExpCoeff * deltaT;

        // Voigt: [σ_xx, σ_yy, σ_zz, σ_xy]
        stress[0] = lambda_d * trace_strain + 2.0 * shearMod * strain_xx - thermalStress - biotCoeff * pressure;
        stress[1] = lambda_d * trace_strain + 2.0 * shearMod * strain_yy - thermalStress - biotCoeff * pressure;
        stress[2] = lambda_d * trace_strain - thermalStress - biotCoeff * pressure; // plane strain: ε_zz = 0
        stress[3] = 2.0 * shearMod * strain_xy;
    } // cauchyStress_PlaneStrain

    // ----------------------------------------------------------------------
    /** Compute fluid content.
     *
     * ζ = α_B ε_v + p/M + 3α_f (T - T_ref)
     */
    static inline
    void fluidContent(const PylithInt dim,
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
                      PylithScalar content[]) {
        const PylithInt i_biotCoeff = numConstants >= 4 ? (PylithInt)constants[0] : 4;
        const PylithInt i_biotMod = numConstants >= 4 ? (PylithInt)constants[1] : 5;
        const PylithInt i_refTemp = numConstants >= 4 ? (PylithInt)constants[2] : 7;
        const PylithInt i_fluidThermalExp = numConstants >= 4 ? (PylithInt)constants[3] : 9;

        const PylithScalar trace_strain = s[sOff[i_trace_strain]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar temperature = s[sOff[i_temperature]];

        const PylithScalar biotCoeff = a[aOff[i_biotCoeff]];
        const PylithScalar biotMod = a[aOff[i_biotMod]];
        const PylithScalar refTemp = a[aOff[i_refTemp]];
        const PylithScalar fluidThermalExp = a[aOff[i_fluidThermalExp]];

        const PylithScalar deltaT = temperature - refTemp;

        content[0] = biotCoeff * trace_strain + pressure / biotMod + 3.0 * fluidThermalExp * deltaT;
    } // fluidContent

    // ----------------------------------------------------------------------
    /** Compute heat flux.
     *
     * q = -k_T ∇T
     */
    static inline
    void heatFlux(const PylithInt dim,
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
                  PylithScalar flux[]) {
        const PylithInt i_thermalCond = numConstants >= 1 ? (PylithInt)constants[0] : 11;

        const PylithScalar* temperature_x = &s_x[sOff_x[i_temperature]];
        const PylithScalar thermalConductivity = a[aOff[i_thermalCond]];

        for (PylithInt d = 0; d < dim; ++d) {
            flux[d] = -thermalConductivity * temperature_x[d];
        }
    } // heatFlux

}; // class IsotropicLinearThermoporoelasticity

// End of file
