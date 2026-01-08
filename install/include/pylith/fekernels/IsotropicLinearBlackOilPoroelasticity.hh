// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity
#include "pylith/fekernels/Tensor.hh" // USES Tensor

#include "pylith/utils/types.hh"

#include <cmath> // USES exp()
#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
/**
 * Kernels for isotropic, linear black oil poroelasticity.
 *
 * The black oil model extends standard poroelasticity with pressure-dependent
 * fluid properties:
 * - Effective viscosity: mu_eff = mu_ref * exp(c_mu * (p - p_ref))
 * - Effective compressibility: c_eff = c_ref * exp(c_c * (p - p_ref))
 */
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticity {
public:

    struct Context {
        PylithReal shearModulus;
        PylithReal drainedBulkModulus;
        PylithReal biotCoefficient;
        PylithReal biotModulus;
        PylithReal referencePressure;
        PylithReal fluidCompressibility;
        PylithReal fluidCompressibilityCoefficient;
        PylithReal viscosityCoefficient;
        Tensor permeability;
        Tensor referenceStress;
        Tensor referenceStrain;
    };

    // --------------------------------------------------------------------------------------------
    static inline
    PylithReal computeEffectiveViscosity(const PylithReal fluidViscosity,
                                         const PylithReal pressure,
                                         const PylithReal referencePressure,
                                         const PylithReal viscosityCoefficient) {
        const PylithReal dp = pressure - referencePressure;
        return fluidViscosity * std::exp(viscosityCoefficient * dp);
    }

    // --------------------------------------------------------------------------------------------
    static inline
    PylithReal computeEffectiveCompressibility(const PylithReal fluidCompressibility,
                                               const PylithReal pressure,
                                               const PylithReal referencePressure,
                                               const PylithReal compressibilityCoefficient) {
        const PylithReal dp = pressure - referencePressure;
        return fluidCompressibility * std::exp(compressibilityCoefficient * dp);
    }

    // --------------------------------------------------------------------------------------------
    static inline
    PylithReal computeEffectiveBiotModulus(const PylithReal porosity,
                                           const PylithReal biotCoefficient,
                                           const PylithReal drainedBulkModulus,
                                           const PylithReal effectiveCompressibility) {
        const PylithReal Ks = (biotCoefficient < 1.0) ?
            drainedBulkModulus / (1.0 - biotCoefficient) : 1.0e50;
        const PylithReal storativity = porosity * effectiveCompressibility +
            (biotCoefficient - porosity) / Ks;
        return (storativity > 0.0) ? 1.0 / storativity : 1.0e50;
    }

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

        // Rheology auxiliary subfields indexed from end of auxiliary array
        // Order: shearModulus, drainedBulkModulus, biotCoefficient, biotModulus,
        //        referencePressure, fluidCompressibility, fluidCompressibilityCoefficient,
        //        viscosityCoefficient, isotropicPermeability
        const PylithInt i_shearModulus = numA - 9;
        const PylithInt i_drainedBulkModulus = numA - 8;
        const PylithInt i_biotCoefficient = numA - 7;
        const PylithInt i_biotModulus = numA - 6;
        const PylithInt i_referencePressure = numA - 5;
        const PylithInt i_fluidCompressibility = numA - 4;
        const PylithInt i_fluidCompressibilityCoefficient = numA - 3;
        const PylithInt i_viscosityCoefficient = numA - 2;
        const PylithInt i_isotropicPermeability = numA - 1;

        assert(numA >= 9);
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_referencePressure] >= 0);
        assert(aOff[i_fluidCompressibility] >= 0);
        assert(aOff[i_fluidCompressibilityCoefficient] >= 0);
        assert(aOff[i_viscosityCoefficient] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);

        context->shearModulus = a[aOff[i_shearModulus]];
        assert(context->shearModulus > 0.0);
        context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        assert(context->drainedBulkModulus > 0.0);
        context->biotCoefficient = a[aOff[i_biotCoefficient]];
        assert(context->biotCoefficient > 0.0);
        context->biotModulus = a[aOff[i_biotModulus]];
        assert(context->biotModulus > 0.0);
        context->referencePressure = a[aOff[i_referencePressure]];
        context->fluidCompressibility = a[aOff[i_fluidCompressibility]];
        assert(context->fluidCompressibility >= 0.0);
        context->fluidCompressibilityCoefficient = a[aOff[i_fluidCompressibilityCoefficient]];
        context->viscosityCoefficient = a[aOff[i_viscosityCoefficient]];

        tensorOps.fromScalar(a[aOff[i_isotropicPermeability]], &context->permeability);
    }

    // --------------------------------------------------------------------------------------------
    static inline
    void cauchyStress(void* rheologyContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(stress);

        const PylithReal shearModulus = context->shearModulus;
        const PylithReal drainedBulkModulus = context->drainedBulkModulus;
        const PylithReal lambda = drainedBulkModulus - 2.0/3.0 * shearModulus;

        const PylithScalar traceStrain = strain.xx + strain.yy + strain.zz;
        const PylithScalar meanStress = lambda * traceStrain;

        stress->xx = meanStress + 2.0 * shearModulus * strain.xx;
        stress->yy = meanStress + 2.0 * shearModulus * strain.yy;
        stress->zz = meanStress + 2.0 * shearModulus * strain.zz;
        stress->xy = 2.0 * shearModulus * strain.xy;
        stress->yz = 2.0 * shearModulus * strain.yz;
        stress->xz = 2.0 * shearModulus * strain.xz;
    }

}; // IsotropicLinearBlackOilPoroelasticity

// ---------------------------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticityPlaneStrain {
public:

    typedef pylith::fekernels::IsotropicLinearBlackOilPoroelasticity BlackOilBase;

    // ================================= LHS =======================================

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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar pressure_t = s_t ? poroContext.pressure_t : 0.0;
        const PylithScalar trace_strain_t = s_t ? poroContext.trace_strain_t : 0.0;

        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, rheologyContext.biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        f0[0] += s_t ? (rheologyContext.biotCoefficient * trace_strain_t) : 0.0;
        f0[0] += s_t ? (pressure_t / effectiveBiotModulus) : 0.0;
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            BlackOilBase::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar* pressure_x = poroContext.pressure_x;
        const PylithScalar fluidViscosity = poroContext.fluidViscosity;

        const PylithReal effectiveViscosity = BlackOilBase::computeEffectiveViscosity(
            fluidViscosity, pressure, rheologyContext.referencePressure, rheologyContext.viscosityCoefficient);

        const PylithReal permeability = rheologyContext.permeability.xx;

        for (PylithInt i = 0; i < dim; ++i) {
            f1[i] += (permeability / effectiveViscosity) * pressure_x[i];
        }
    }

    // ================================= JACOBIANS =======================================

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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithReal shearModulus = rheologyContext.shearModulus;

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            }
        }
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf2[i * _dim + i] += -biotCoefficient;
        }
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithReal drainedBulkModulus = rheologyContext.drainedBulkModulus;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf2[i * _dim + i] += drainedBulkModulus;
        }
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, rheologyContext.biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        Jf0[0] += s_tshift / effectiveBiotModulus;
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar fluidViscosity = poroContext.fluidViscosity;

        const PylithReal effectiveViscosity = BlackOilBase::computeEffectiveViscosity(
            fluidViscosity, pressure, rheologyContext.referencePressure, rheologyContext.viscosityCoefficient);

        const PylithReal permeability = rheologyContext.permeability.xx;
        const PylithReal darcyConductivity = permeability / effectiveViscosity;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf3[i * _dim + i] += darcyConductivity;
        }
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        Jf0[0] += rheologyContext.biotCoefficient * s_tshift;
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        // First add the standard f0p_implicit terms
        f0p_implicit(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                     t, x, numConstants, constants, f0);

        // Add source density term from Poroelasticity context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        f0[0] += poroContext.sourceDensity;
    }

    // ================================= DERIVED FIELDS =======================================

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
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain,
            BlackOilBase::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            stressVector);
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar traceStrain = poroContext.trace_strain;
        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;

        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        waterContent[0] = biotCoefficient * traceStrain + pressure / effectiveBiotModulus;
    }

}; // IsotropicLinearBlackOilPoroelasticityPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
class pylith::fekernels::IsotropicLinearBlackOilPoroelasticity3D {
public:

    typedef pylith::fekernels::IsotropicLinearBlackOilPoroelasticity BlackOilBase;

    // ================================= LHS =======================================

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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar pressure_t = s_t ? poroContext.pressure_t : 0.0;
        const PylithScalar trace_strain_t = s_t ? poroContext.trace_strain_t : 0.0;

        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, rheologyContext.biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        f0[0] += s_t ? (rheologyContext.biotCoefficient * trace_strain_t) : 0.0;
        f0[0] += s_t ? (pressure_t / effectiveBiotModulus) : 0.0;
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::f1v(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            BlackOilBase::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            f1);
    }

    // --------------------------------------------------------------------------------------------
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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar* pressure_x = poroContext.pressure_x;
        const PylithScalar fluidViscosity = poroContext.fluidViscosity;

        const PylithReal effectiveViscosity = BlackOilBase::computeEffectiveViscosity(
            fluidViscosity, pressure, rheologyContext.referencePressure, rheologyContext.viscosityCoefficient);

        const PylithReal permeability = rheologyContext.permeability.xx;

        for (PylithInt i = 0; i < dim; ++i) {
            f1[i] += (permeability / effectiveViscosity) * pressure_x[i];
        }
    }

    // ================================= JACOBIANS =======================================

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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithReal shearModulus = rheologyContext.shearModulus;
        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            }
        }
    }

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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf2[i * _dim + i] += -biotCoefficient;
        }
    }

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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithReal drainedBulkModulus = rheologyContext.drainedBulkModulus;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf2[i * _dim + i] += drainedBulkModulus;
        }
    }

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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, rheologyContext.biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        Jf0[0] += s_tshift / effectiveBiotModulus;
    }

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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar fluidViscosity = poroContext.fluidViscosity;

        const PylithReal effectiveViscosity = BlackOilBase::computeEffectiveViscosity(
            fluidViscosity, pressure, rheologyContext.referencePressure, rheologyContext.viscosityCoefficient);

        const PylithReal permeability = rheologyContext.permeability.xx;
        const PylithReal darcyConductivity = permeability / effectiveViscosity;

        for (PylithInt i = 0; i < _dim; ++i) {
            Jf3[i * _dim + i] += darcyConductivity;
        }
    }

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
        assert(_dim == dim);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        Jf0[0] += rheologyContext.biotCoefficient * s_tshift;
    }

    // --------------------------------------------------------------------------------------------
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
        const PylithInt _dim = 3;
        assert(_dim == dim);

        // First add the standard f0p_implicit terms
        f0p_implicit(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                     t, x, numConstants, constants, f0);

        // Add source density term from Poroelasticity context
        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextSourceDensity(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        f0[0] += poroContext.sourceDensity;
    }

    // ================================= DERIVED FIELDS =======================================

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
        assert(_dim == dim);

        pylith::fekernels::Elasticity::StrainContext strainContext;
        pylith::fekernels::Elasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Elasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Elasticity3D::infinitesimalStrain,
            BlackOilBase::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            stressVector);
    }

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
        assert(_dim == dim);

        pylith::fekernels::Poroelasticity::Context poroContext;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &poroContext, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        BlackOilBase::Context rheologyContext;
        BlackOilBase::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        const PylithScalar pressure = poroContext.pressure;
        const PylithScalar traceStrain = poroContext.trace_strain;
        const PylithReal biotCoefficient = rheologyContext.biotCoefficient;

        const PylithReal effectiveCompressibility = BlackOilBase::computeEffectiveCompressibility(
            rheologyContext.fluidCompressibility, pressure,
            rheologyContext.referencePressure, rheologyContext.fluidCompressibilityCoefficient);

        const PylithReal effectiveBiotModulus = BlackOilBase::computeEffectiveBiotModulus(
            poroContext.porosity, biotCoefficient,
            rheologyContext.drainedBulkModulus, effectiveCompressibility);

        waterContent[0] = biotCoefficient * traceStrain + pressure / effectiveBiotModulus;
    }

}; // IsotropicLinearBlackOilPoroelasticity3D

// End of file
