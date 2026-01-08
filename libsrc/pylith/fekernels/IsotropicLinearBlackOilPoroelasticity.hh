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
 * Kernels for isotropic linear poroelasticity with black oil formulation.
 *
 * The black oil model extends standard poroelasticity with pressure-dependent
 * fluid properties. The key extensions include:
 *
 * 1. Pressure-dependent fluid compressibility:
 *    c_f(p) = c_f0 + c_f1 * (p - p_ref)
 *    where c_f0 is reference compressibility, c_f1 is the pressure coefficient
 *
 * 2. Pressure-dependent fluid viscosity:
 *    mu(p) = mu_ref * exp(alpha_mu * (p - p_ref))
 *    where mu_ref is reference viscosity, alpha_mu is the viscosity coefficient
 *
 * Solution fields: [disp(dim), pressure(1), trace_strain(1)] (QS)
 *
 * Auxiliary fields:
 * - 0: solid_density(1)
 * - 1: fluid_density(1)
 * - 2: fluid_viscosity(1) - reference viscosity
 * - 3: porosity(1)
 *
 * Optional fields:
 * - +1: gravity_field (dim, optional)
 * - +1: body_force(dim, optional)
 * - +1: source_density(1, optional)
 * - +1: reference_stress (optional)
 * - +1: reference_strain (optional)
 *
 * Rheology fields (black oil specific):
 * - numA - 9: shear_modulus(1)
 * - numA - 8: drained_bulk_modulus(1)
 * - numA - 7: biot_coefficient(1)
 * - numA - 6: biot_modulus(1)
 * - numA - 5: reference_pressure(1)
 * - numA - 4: fluid_compressibility(1)
 * - numA - 3: fluid_compressibility_coefficient(1)
 * - numA - 2: viscosity_coefficient(1)
 * - numA - 1: isotropic_permeability(1) OR tensor_permeability
 */

#include "pylith/fekernels/fekernelsfwd.hh"
#include "pylith/fekernels/Poroelasticity.hh"
#include "pylith/fekernels/Elasticity.hh"

#include "pylith/utils/types.hh"

#include <cassert>

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity with black oil formulation.
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticity {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal pressure;
        PylithReal trace_strain;
        PylithReal trace_strain_t;
        PylithReal fluidViscosity;        // Reference viscosity
        PylithReal shearModulus;
        PylithReal drainedBulkModulus;
        PylithReal biotCoefficient;
        PylithReal biotModulus;
        PylithReal referencePressure;      // Black oil: reference pressure
        PylithReal fluidCompressibility;   // Black oil: reference compressibility
        PylithReal fluidCompressibilityCoefficient; // Black oil: pressure dependence of compressibility
        PylithReal viscosityCoefficient;   // Black oil: pressure dependence of viscosity
        pylith::fekernels::Tensor permeability;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Compute pressure-dependent effective viscosity for black oil model.
     *
     * mu_eff(p) = mu_ref * exp(alpha_mu * (p - p_ref))
     */
    static inline
    PylithReal computeEffectiveViscosity(const Context& context) {
        const PylithReal deltaP = context.pressure - context.referencePressure;
        // Use linear approximation for small deltaP, full exponential otherwise
        const PylithReal factor = context.viscosityCoefficient * deltaP;
        if (fabs(factor) < 0.01) {
            return context.fluidViscosity * (1.0 + factor);
        }
        return context.fluidViscosity * exp(factor);
    } // computeEffectiveViscosity

    // --------------------------------------------------------------------------------------------
    /** Compute pressure-dependent effective compressibility for black oil model.
     *
     * c_f(p) = c_f0 + c_f1 * (p - p_ref)
     */
    static inline
    PylithReal computeEffectiveCompressibility(const Context& context) {
        const PylithReal deltaP = context.pressure - context.referencePressure;
        PylithReal c_eff = context.fluidCompressibility +
                          context.fluidCompressibilityCoefficient * deltaP;
        // Ensure non-negative compressibility
        return (c_eff > 0.0) ? c_eff : context.fluidCompressibility;
    } // computeEffectiveCompressibility

    // --------------------------------------------------------------------------------------------
    /** Compute effective Biot modulus with pressure-dependent fluid compressibility.
     *
     * 1/M_eff = 1/M_solid + phi * c_f(p)
     *
     * where M_solid is the original Biot modulus (without fluid compressibility)
     * and c_f(p) is the pressure-dependent fluid compressibility.
     */
    static inline
    PylithReal computeEffectiveBiotModulus(const Context& context,
                                           const PylithReal porosity) {
        const PylithReal c_f = computeEffectiveCompressibility(context);
        const PylithReal oneOverM_solid = 1.0 / context.biotModulus;
        const PylithReal oneOverM_eff = oneOverM_solid + porosity * c_f;
        return 1.0 / oneOverM_eff;
    } // computeEffectiveBiotModulus

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext(Context* context,
                    const PylithInt dim,
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
                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming solution subfields
        const PylithInt i_pressure = 1;

        // Incoming poroelastic auxiliary subfields
        const PylithInt i_fluidViscosity = 2;

        // Incoming rheology auxiliary subfields (black oil)
        const PylithInt i_shearModulus = numA - 9;
        const PylithInt i_drainedBulkModulus = numA - 8;
        const PylithInt i_biotCoefficient = numA - 7;
        const PylithInt i_biotModulus = numA - 6;
        const PylithInt i_referencePressure = numA - 5;
        const PylithInt i_fluidCompressibility = numA - 4;
        const PylithInt i_fluidCompressibilityCoefficient = numA - 3;
        const PylithInt i_viscosityCoefficient = numA - 2;

        assert(numA >= 9);
        assert(s);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_referencePressure] >= 0);
        assert(aOff[i_fluidCompressibility] >= 0);
        assert(aOff[i_fluidCompressibilityCoefficient] >= 0);
        assert(aOff[i_viscosityCoefficient] >= 0);

        // Solution Variables
        context->pressure = s[sOff[i_pressure]];

        // Poroelastic Auxiliary Variables
        context->fluidViscosity = a[aOff[i_fluidViscosity]];
        assert(context->fluidViscosity > 0.0);

        // Rheology Specific Auxiliary Variables
        context->shearModulus = a[aOff[i_shearModulus]];
        assert(context->shearModulus > 0.0);
        context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        assert(context->drainedBulkModulus > 0.0);
        context->biotCoefficient = a[aOff[i_biotCoefficient]];
        assert(context->biotCoefficient > 0.0);
        context->biotModulus = a[aOff[i_biotModulus]];
        assert(context->biotModulus > 0.0);

        // Black Oil Specific Variables
        context->referencePressure = a[aOff[i_referencePressure]];
        context->fluidCompressibility = a[aOff[i_fluidCompressibility]];
        context->fluidCompressibilityCoefficient = a[aOff[i_fluidCompressibilityCoefficient]];
        context->viscosityCoefficient = a[aOff[i_viscosityCoefficient]];

    } // setContext

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextIsotropicPerm(Context* context,
                                 const PylithInt dim,
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
                                 const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming auxiliary fields.
        const PylithInt i_isotropicPermeability = numA - 1;

        // Using isotropic permeability
        tensorOps.fromScalar(a[aOff[i_isotropicPermeability]], &context->permeability);

    } // setContextIsotropicPerm

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextTensorPerm(Context* context,
                              const PylithInt dim,
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
                              const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        // Incoming auxiliary fields.
        const PylithInt i_tensorPermeability = numA - 1;

        // Using tensor permeability
        tensorOps.fromVector(&a[aOff[i_tensorPermeability]], &context->permeability);

    } // setContextTensorPerm

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextQuasistatic(Context* context,
                               const PylithInt dim,
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
                               const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);
        // Incoming solution fields.
        const PylithInt i_trace_strain = 2;
        assert(sOff[i_trace_strain] >= 0);

        // Variables &c
        context->trace_strain = s[sOff[i_trace_strain]];

    } // setContextQuasistatic

    // --------------------------------------------------------------------------------------------
    static inline
    void setContextRefState(Context* context,
                            const PylithInt dim,
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
                            const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);
        // Incoming auxiliary fields (with reference state, they come before rheology fields)
        const PylithInt i_refStress = numA - 11;
        const PylithInt i_refStrain = numA - 10;

        // Reference stress and strain
        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);

    } // setContextRefState

    // --------------------------------------------------------------------------------------------
    /** Compute effective Darcy conductivity with pressure-dependent viscosity.
     *
     * K = k / mu_eff(p)
     */
    static inline
    void darcyConductivityBlackOil(const Context& context,
                                   const TensorOps& tensorOps,
                                   Tensor* darcyConductivity) {
        assert(darcyConductivity);

        const PylithReal mu_eff = computeEffectiveViscosity(context);
        const PylithReal conductivityScale = 1.0 / mu_eff;

        tensorOps.scale(conductivityScale, context.permeability, darcyConductivity);
    } // darcyConductivityBlackOil

    // --------------------------------------------------------------------------------------------
    /** Compute mean stress from trace strain and pressure for black oil model.
     */
    static inline
    void meanStress(const Context& context,
                    pylith::fekernels::Scalar* meanStress) {
        assert(meanStress);

        const PylithReal lambda = context.drainedBulkModulus - 2.0/3.0*context.shearModulus;
        const PylithReal K_d = context.drainedBulkModulus;

        meanStress->value = K_d * context.trace_strain - context.biotCoefficient * context.pressure;
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Compute deviatoric stress from strain for black oil model.
     */
    static inline
    void deviatoricStress(const Context& context,
                          const Tensor& strain,
                          const TensorOps& tensorOps,
                          Tensor* deviatoricStress) {
        assert(deviatoricStress);

        tensorOps.deviator(strain, deviatoricStress);
        tensorOps.scale(2.0 * context.shearModulus, *deviatoricStress, deviatoricStress);
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Compute Cauchy stress from strain for black oil model.
     */
    static inline
    void cauchyStress(const Context& context,
                      const Tensor& strain,
                      const TensorOps& tensorOps,
                      Tensor* cauchyStress) {
        assert(cauchyStress);

        Scalar meanStressScalar;
        meanStress(context, &meanStressScalar);

        Tensor deviatoricStressTensor;
        deviatoricStress(context, strain, tensorOps, &deviatoricStressTensor);

        tensorOps.addScalarToTrace(meanStressScalar, deviatoricStressTensor, cauchyStress);
    } // cauchyStress

}; // IsotropicLinearBlackOilPoroelasticity


// ================================================================================================
// Kernels for isotropic, linear poroelasticity with black oil formulation - Plane strain
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // Use 2D tensor operations
    typedef pylith::fekernels::IsotropicLinearBlackOilPoroelasticity BlackOilBase;

    // ============================= LHS Residual ================================================

    // --------------------------------------------------------------------------------------------
    /** f0p function for implicit time stepping (fluid mass balance equation).
     *
     * f0p = \dot{zeta} - source
     *     = (1/M_eff) * dp/dt + alpha * d(epsilon_v)/dt - source
     *
     * With black oil: M_eff includes pressure-dependent fluid compressibility
     */
    static inline
    void f0p_implicit(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute effective Biot modulus with pressure-dependent compressibility
        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);
        const PylithReal alpha = rheologyContext.biotCoefficient;

        // Time derivatives
        const PylithReal pressure_t = poroContext.pressure_t;
        const PylithReal trace_strain_t = rheologyContext.trace_strain_t;

        // f0p = (1/M_eff) * dp/dt + alpha * d(epsilon_v)/dt
        f0[0] += pressure_t / M_eff + alpha * trace_strain_t;

    } // f0p_implicit

    // --------------------------------------------------------------------------------------------
    /** f0p function for implicit time stepping with source density.
     */
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Call base function
        f0p_implicit(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                     aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0);

        // Add source density
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        f0[0] -= poroContext.sourceDensity;

    } // f0p_implicit_source

    // --------------------------------------------------------------------------------------------
    /** f1u function for implicit time stepping (momentum equation).
     *
     * f1u = stress tensor
     */
    static inline
    void f1u(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute strain from displacement gradient
        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity::infinitesimalStrain(poroContext.displacement_x, tensorOps, &strain);

        // Compute Cauchy stress
        pylith::fekernels::Tensor stress;
        BlackOilBase::cauchyStress(rheologyContext, strain, tensorOps, &stress);

        // Convert to f1
        PylithReal stressTensor[4] = {0.0, 0.0, 0.0, 0.0};
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] += stressTensor[i];
        }

    } // f1u

    // --------------------------------------------------------------------------------------------
    /** f1p function for implicit time stepping (Darcy flux).
     *
     * f1p = -K * grad(p) with pressure-dependent viscosity
     */
    static inline
    void f1p(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextIsotropicPerm(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute Darcy conductivity with pressure-dependent viscosity
        pylith::fekernels::Tensor darcyConductivity;
        BlackOilBase::darcyConductivityBlackOil(rheologyContext, tensorOps, &darcyConductivity);

        // Compute flux: -K * grad(p)
        const PylithScalar* pressure_x = poroContext.pressure_x;
        PylithReal conductivityTensor[4] = {0.0, 0.0, 0.0, 0.0};
        tensorOps.toTensor(darcyConductivity, conductivityTensor);

        for (PylithInt i = 0; i < _dim; ++i) {
            f1[i] = 0.0;
            for (PylithInt j = 0; j < _dim; ++j) {
                f1[i] -= conductivityTensor[i*_dim+j] * pressure_x[j];
            }
        }

    } // f1p

    // ============================= LHS Jacobian ================================================

    // --------------------------------------------------------------------------------------------
    /** Jf3uu function - elastic constants contribution to Jacobian.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal G = rheologyContext.shearModulus;
        const PylithReal K_d = rheologyContext.drainedBulkModulus;
        const PylithReal lambda = K_d - 2.0/3.0*G;

        // Same as standard poroelasticity
        const PylithReal C1111 = lambda + 2.0*G;
        const PylithReal C1122 = lambda;
        const PylithReal C1212 = G;

        /* j(f,g,df,dg) = C(f,df,g,dg)
         *
         * 0:  j0000 = C1111 = lambda + 2*G
         * 1:  j0001 = C1112 = 0
         * 2:  j0010 = C1121 = 0
         * 3:  j0011 = C1122 = lambda
         * 4:  j0100 = C1211 = 0
         * 5:  j0101 = C1212 = G
         * 6:  j0110 = C1221 = G
         * 7:  j0111 = C1222 = 0
         * 8:  j1000 = C2111 = 0
         * 9:  j1001 = C2112 = G
         * 10: j1010 = C2121 = G
         * 11: j1011 = C2122 = 0
         * 12: j1100 = C2211 = lambda
         * 13: j1101 = C2212 = 0
         * 14: j1110 = C2221 = 0
         * 15: j1111 = C2222 = lambda + 2*G
         */
        Jf3[ 0] -= C1111; // j0000
        Jf3[ 3] -= C1122; // j0011
        Jf3[ 5] -= C1212; // j0101
        Jf3[ 6] -= C1212; // j0110
        Jf3[ 9] -= C1212; // j1001
        Jf3[10] -= C1212; // j1010
        Jf3[12] -= C1122; // j1100
        Jf3[15] -= C1111; // j1111

    } // Jf3uu

    // --------------------------------------------------------------------------------------------
    /** Jf2up function - Biot coefficient contribution to Jacobian.
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
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal alpha = rheologyContext.biotCoefficient;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d*_dim+d] += alpha;
        }

    } // Jf2up

    // --------------------------------------------------------------------------------------------
    /** Jf2ue function - Lambda contribution to Jacobian.
     */
    static inline
    void Jf2ue(const PylithInt dim,
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
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal K_d = rheologyContext.drainedBulkModulus;
        const PylithReal G = rheologyContext.shearModulus;
        const PylithReal lambda = K_d - 2.0/3.0*G;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d*_dim+d] -= lambda;
        }

    } // Jf2ue

    // --------------------------------------------------------------------------------------------
    /** Jf0pp function - Storage coefficient contribution to Jacobian.
     *
     * For black oil model, this includes pressure-dependent compressibility.
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
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context for porosity
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute effective Biot modulus with pressure-dependent compressibility
        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);

        Jf0[0] += s_tshift / M_eff;

    } // Jf0pp

    // --------------------------------------------------------------------------------------------
    /** Jf3pp function - Darcy conductivity contribution to Jacobian.
     *
     * For black oil model, includes pressure-dependent viscosity.
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
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextIsotropicPerm(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute Darcy conductivity with pressure-dependent viscosity
        pylith::fekernels::Tensor darcyConductivity;
        BlackOilBase::darcyConductivityBlackOil(rheologyContext, tensorOps, &darcyConductivity);

        PylithReal conductivityTensor[4] = {0.0, 0.0, 0.0, 0.0};
        tensorOps.toTensor(darcyConductivity, conductivityTensor);

        // Jf3[i*dim*dim*dim + j*dim*dim + k*dim + l] corresponds to dF_i/d(grad_l p_k) * dgrad_j
        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[i*_dim + j] -= conductivityTensor[i*_dim+j];
            }
        }

    } // Jf3pp

    // --------------------------------------------------------------------------------------------
    /** Jf0pe function - Biot coefficient contribution to pressure-strain coupling.
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
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal alpha = rheologyContext.biotCoefficient;

        Jf0[0] += s_tshift * alpha;

    } // Jf0pe

    // ============================= Derived Fields ==============================================

    // --------------------------------------------------------------------------------------------
    /** Compute Cauchy stress for derived field.
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute strain from displacement gradient
        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity::infinitesimalStrain(poroContext.displacement_x, tensorOps, &strain);

        // Compute Cauchy stress
        pylith::fekernels::Tensor stress;
        BlackOilBase::cauchyStress(rheologyContext, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Compute water content for derived field.
     */
    static inline
    void waterContent_asScalar(const PylithInt dim,
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
                               PylithScalar waterContent[]) {
        const PylithInt _dim = 2;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute effective Biot modulus with pressure-dependent compressibility
        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);
        const PylithReal alpha = rheologyContext.biotCoefficient;

        // Water content = p/M_eff + alpha * epsilon_v
        waterContent[0] = rheologyContext.pressure / M_eff + alpha * rheologyContext.trace_strain;

    } // waterContent_asScalar

}; // IsotropicLinearBlackOilPoroelasticityPlaneStrain


// ================================================================================================
// Kernels for isotropic, linear poroelasticity with black oil formulation - 3D
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    typedef pylith::fekernels::IsotropicLinearBlackOilPoroelasticity BlackOilBase;

    // ============================= LHS Residual ================================================

    // --------------------------------------------------------------------------------------------
    /** f0p function for implicit time stepping.
     */
    static inline
    void f0p_implicit(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Black oil rheology context
        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        // Compute effective Biot modulus with pressure-dependent compressibility
        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);
        const PylithReal alpha = rheologyContext.biotCoefficient;

        // Time derivatives
        const PylithReal pressure_t = poroContext.pressure_t;
        const PylithReal trace_strain_t = rheologyContext.trace_strain_t;

        f0[0] += pressure_t / M_eff + alpha * trace_strain_t;

    } // f0p_implicit

    // --------------------------------------------------------------------------------------------
    /** f0p function with source.
     */
    static inline
    void f0p_implicit_source(const PylithInt dim,
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
        f0p_implicit(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                     aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, f0);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        f0[0] -= poroContext.sourceDensity;

    } // f0p_implicit_source

    // --------------------------------------------------------------------------------------------
    /** f1u function for 3D.
     */
    static inline
    void f1u(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity::infinitesimalStrain(poroContext.displacement_x, tensorOps, &strain);

        pylith::fekernels::Tensor stress;
        BlackOilBase::cauchyStress(rheologyContext, strain, tensorOps, &stress);

        PylithReal stressTensor[9] = {0.0};
        tensorOps.toTensor(stress, stressTensor);

        for (PylithInt i = 0; i < _dim*_dim; ++i) {
            f1[i] += stressTensor[i];
        }

    } // f1u

    // --------------------------------------------------------------------------------------------
    /** f1p function for 3D.
     */
    static inline
    void f1p(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextIsotropicPerm(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        pylith::fekernels::Tensor darcyConductivity;
        BlackOilBase::darcyConductivityBlackOil(rheologyContext, tensorOps, &darcyConductivity);

        const PylithScalar* pressure_x = poroContext.pressure_x;
        PylithReal conductivityTensor[9] = {0.0};
        tensorOps.toTensor(darcyConductivity, conductivityTensor);

        for (PylithInt i = 0; i < _dim; ++i) {
            f1[i] = 0.0;
            for (PylithInt j = 0; j < _dim; ++j) {
                f1[i] -= conductivityTensor[i*_dim+j] * pressure_x[j];
            }
        }

    } // f1p

    // ============================= LHS Jacobian ================================================

    // --------------------------------------------------------------------------------------------
    /** Jf3uu function for 3D.
     */
    static inline
    void Jf3uu(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal G = rheologyContext.shearModulus;
        const PylithReal K_d = rheologyContext.drainedBulkModulus;
        const PylithReal lambda = K_d - 2.0/3.0*G;

        const PylithReal C1111 = lambda + 2.0*G;
        const PylithReal C1122 = lambda;
        const PylithReal C1212 = G;

        // 3D elasticity tensor - standard poroelasticity form
        /* Nonzero entries in order:
         * C1111, C1122, C1133, C1212, C1313, C2222, C2233, C2323, C3333
         */
        Jf3[ 0] -= C1111; // C1111
        Jf3[ 4] -= C1122; // C1122
        Jf3[ 8] -= C1122; // C1133
        Jf3[10] -= C1212; // C1212
        Jf3[12] -= C1212; // C1221
        Jf3[20] -= C1212; // C1313
        Jf3[24] -= C1212; // C1331
        Jf3[28] -= C1212; // C2112
        Jf3[30] -= C1212; // C2121
        Jf3[36] -= C1122; // C2211
        Jf3[40] -= C1111; // C2222
        Jf3[44] -= C1122; // C2233
        Jf3[50] -= C1212; // C2323
        Jf3[52] -= C1212; // C2332
        Jf3[56] -= C1212; // C3113
        Jf3[60] -= C1212; // C3131
        Jf3[68] -= C1212; // C3223
        Jf3[70] -= C1212; // C3232
        Jf3[72] -= C1122; // C3311
        Jf3[76] -= C1122; // C3322
        Jf3[80] -= C1111; // C3333

    } // Jf3uu

    // --------------------------------------------------------------------------------------------
    /** Jf2up function for 3D.
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
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal alpha = rheologyContext.biotCoefficient;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d*_dim+d] += alpha;
        }

    } // Jf2up

    // --------------------------------------------------------------------------------------------
    /** Jf2ue function for 3D.
     */
    static inline
    void Jf2ue(const PylithInt dim,
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
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal K_d = rheologyContext.drainedBulkModulus;
        const PylithReal G = rheologyContext.shearModulus;
        const PylithReal lambda = K_d - 2.0/3.0*G;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d*_dim+d] -= lambda;
        }

    } // Jf2ue

    // --------------------------------------------------------------------------------------------
    /** Jf0pp function for 3D.
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
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);

        Jf0[0] += s_tshift / M_eff;

    } // Jf0pp

    // --------------------------------------------------------------------------------------------
    /** Jf3pp function for 3D.
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
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextIsotropicPerm(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                              aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        pylith::fekernels::Tensor darcyConductivity;
        BlackOilBase::darcyConductivityBlackOil(rheologyContext, tensorOps, &darcyConductivity);

        PylithReal conductivityTensor[9] = {0.0};
        tensorOps.toTensor(darcyConductivity, conductivityTensor);

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[i*_dim + j] -= conductivityTensor[i*_dim+j];
            }
        }

    } // Jf3pp

    // --------------------------------------------------------------------------------------------
    /** Jf0pe function for 3D.
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
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal alpha = rheologyContext.biotCoefficient;

        Jf0[0] += s_tshift * alpha;

    } // Jf0pe

    // ============================= Derived Fields ==============================================

    // --------------------------------------------------------------------------------------------
    /** Compute Cauchy stress for derived field.
     */
    static inline
    void cauchyStress_infinitesimalStrain_asVector(const PylithInt dim,
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
        const pylith::fekernels::TensorOps tensorOps(_dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        pylith::fekernels::Tensor strain;
        pylith::fekernels::Elasticity::infinitesimalStrain(poroContext.displacement_x, tensorOps, &strain);

        pylith::fekernels::Tensor stress;
        BlackOilBase::cauchyStress(rheologyContext, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Compute water content for derived field.
     */
    static inline
    void waterContent_asScalar(const PylithInt dim,
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
                               PylithScalar waterContent[]) {
        const PylithInt _dim = 3;
        const pylith::fekernels::TensorOps tensorOps(_dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                 aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);
        BlackOilBase::setContextQuasistatic(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, tensorOps);

        const PylithReal M_eff = BlackOilBase::computeEffectiveBiotModulus(rheologyContext, poroContext.porosity);
        const PylithReal alpha = rheologyContext.biotCoefficient;

        waterContent[0] = rheologyContext.pressure / M_eff + alpha * rheologyContext.trace_strain;

    } // waterContent_asScalar

}; // IsotropicLinearBlackOilPoroelasticity3D


// End of file
