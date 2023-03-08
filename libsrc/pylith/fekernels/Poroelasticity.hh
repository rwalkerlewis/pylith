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
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/Poroelasticity.hh
 *
 * Solution fields: [disp(dim), pres, vol_strain]
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * /** Kernel interface.
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
 * @param[out] Storage coefficient at constant strain.
 */

#if !defined(pylith_fekernels_poroelasticity_hh)
#define pylith_fekernels_poroelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Tensor.hh"

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::Poroelasticity {
public:

    struct StrainContext {
        PylithInt dim;
        const PylithReal* disp;
        const PylithReal* disp_t;
        const PylithReal* disp_x;
        const PylithReal* x;
    };

    // Interface for functions computing strain.
    typedef void (*strainfn_type) (const StrainContext& context,
                                   pylith::fekernels::Tensor*);

    // Interface for functions computing stress.
    typedef void (*stressfn_type) (void*,
                                   const pylith::fekernels::Tensor&,
                                   const pylith::fekernels::TensorOps&,
                                   pylith::fekernels::Tensor*);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Set strain context.
     */
    static inline
    void setStrainContext(StrainContext* context,
                          const PylithInt dim,
                          const PylithInt numS,
                          const PylithInt sOff[],
                          const PylithInt sOff_x[],
                          const PylithScalar s[],
                          const PylithScalar s_t[],
                          const PylithScalar s_x[],
                          const PylithScalar x[]) {
        assert(context);
        assert(numS >= 1);

        const PylithInt i_disp = 0;

        assert(sOff[i_disp] >= 0);

        context->dim = dim;
        context->disp = &s[sOff[i_disp]];
        context->disp_t = &s_t[sOff[i_disp]];
        context->disp_x = &s_x[sOff_x[i_disp]];
        context->x = x;
    }

    // =============================================================================
    // Displacement
    // =============================================================================
    // ----------------------------------------------------------------------
    /** f0 function for poroelasticity equation for displacement field
     *
     * Placeholder function
     *
     */
    static inline
    void f0u(const PylithInt dim,
             const PylithInt numS,
             const PylithInt numA,
             const PylithInt sOff[],
             const PylithInt sOff_x[],
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
        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += 0.0;
        } // for
    } // f0u

    // =============================================================================
    // Velocity
    // =============================================================================
    // ----------------------------------------------------------------------
    /** f0 function for implicit time stepping pporoelasticity equation for velocity field
     *
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
        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_velocity = 3;

        assert(sOff);
        assert(sOff[i_displacement] >= 0);
        assert(sOff[i_velocity] >= 0);
        assert(s_t);
        assert(aOff);
        assert(s);

        const PylithScalar *displacement_t = &s_t[sOff[i_displacement]]; // disp_t
        const PylithScalar *velocity = &s[sOff[i_velocity]]; // vel

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += displacement_t[i];
            f0[i] -= velocity[i];
        } // for
    } // f0v_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * f0v function for poroelasticity equation, explicit time stepping, dynamic.
     *
     */

    static inline
    void f0v_explicit(const PylithInt dim,
                      const PylithInt numS,
                      const PylithInt numA,
                      const PylithInt sOff[],
                      const PylithInt sOff_x[],
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
        // Incoming solution fields.
        const PylithInt i_velocity = 2;

        // Incoming auxiliary fields.
        const PylithInt i_solid_density = 0;
        const PylithInt i_fluid_density = 1;
        const PylithInt i_porosity = 3;

        const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
        const PylithScalar* velocity_t = &s_t[sOff[i_velocity]]; // acceleration

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += velocity_t[i] * bulkDensity;
        } // for
    } // f0v_explicit

    // =============================================================================
    // Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    // f0e function for isotropic linear Poroelasticity.
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
        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_trace_strain = 2;

        // Incoming auxiliary fields.

        const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        for (PylithInt d = 0; d < dim; ++d) {
            f0[0] += displacement_x[d*dim+d];
        }
        f0[0] -= trace_strain;
    } // f0e

    // =============================================================================
    // Time Derivative of Pressure
    // =============================================================================
    // ----------------------------------------------------------------------
    /*
     * f0 function for poroelasticity equation, implicit time stepping, quasistatic, for
     * time derivative of pressure.
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
        // Incoming solution fields.
        const PylithInt i_pressure = 1;
        const PylithInt i_pdot = 4;

        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_pdot] >= 0);
        assert(s);
        assert(s_t);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]]; // disp_t
        const PylithScalar pdot = s[sOff[i_pdot]]; // vel

        f0[0] += pressure_t;
        f0[0] -= pdot;
    } // f0pdot

    // =============================================================================
    // Time Derivative of Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    /*
     * f0 function for poroelasticity equation, implicit time stepping, quasistatic, for
     * time derivative of volumetric strain.
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
        // Incoming solution fields.
        const PylithInt i_trace_strain = 2;
        const PylithInt i_edot = 5;

        assert(sOff);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff[i_edot] >= 0);

        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]]; // disp_t
        const PylithScalar edot = s[sOff[i_edot]]; // vel

        f0[0] += trace_strain_t - edot;
    } // f0edot

    /* -------------------------------------------------------------------------- */
    /*                           RHS Residuals                                    */
    /* -------------------------------------------------------------------------- */
    // Quasi-Static

    // =============================================================================
    // Displacement
    // =============================================================================

    // ----------------------------------------------------------------------
    /*
     * g0 function for displacement equation: g0u = v.
     */
    static inline
    void g0u(const PylithInt dim,
             const PylithInt numS,
             const PylithInt numA,
             const PylithInt sOff[],
             const PylithInt sOff_x[],
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
        const PylithInt i_velocity = 2;
        const PylithScalar* velocity = &s[sOff[i_velocity]];

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += velocity[i];
        } // for
    } // g0u

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
     */
    static inline
    void g0v_grav(const PylithInt dim,
                  const PylithInt numS,
                  const PylithInt numA,
                  const PylithInt sOff[],
                  const PylithInt sOff_x[],
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
        // Incoming auxililary fields.

        // Poroelasticity
        const PylithInt i_solid_density = 0;
        const PylithInt i_fluid_density = 1;
        const PylithInt i_porosity = 3;

        // 3 + n
        const PylithInt i_gravityField = 4;

        const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
        const PylithScalar* gravityField = &a[aOff[i_gravityField]];

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bulkDensity * gravityField[i];
        } // for
    } // g0v_grav

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
     */
    static inline
    void g0v_bodyforce(const PylithInt dim,
                       const PylithInt numS,
                       const PylithInt numA,
                       const PylithInt sOff[],
                       const PylithInt sOff_x[],
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
        // Incoming auxiliary fields

        // Poroelasticity

        // 3 + n
        const PylithInt i_bodyForce = 4;

        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bodyForce[i];
        } // for
    } // g0v_bodyforce

    // ----------------------------------------------------------------------
    /*
     * g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
     */
    static inline
    void g0v_grav_bodyforce(const PylithInt dim,
                            const PylithInt numS,
                            const PylithInt numA,
                            const PylithInt sOff[],
                            const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity
        const PylithInt i_solid_density = 0;
        const PylithInt i_fluid_density = 1;
        const PylithInt i_porosity = 3;

        // 3 + n
        const PylithInt i_bodyForce = 4;
        const PylithInt i_gravityField = 5;

        const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
        const PylithScalar* gravityField = &a[aOff[i_gravityField]];
        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

        // gravity field
        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bulkDensity * gravityField[i];
        } // for

        // body force
        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bodyForce[i];
        } // for
    } // g0v_grav_bodyforce

    // =============================================================================
    // Pressure
    // =============================================================================

    // ----------------------------------------------------------------------
    /* g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
     */
    static inline
    void g0p_source(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity

        const PylithInt i_sourceDensity = 0;
        const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += sourceDensity[i];
        } // for
    } // g0p_source

    // ----------------------------------------------------------------------
    /* g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
     *
     */
    static inline
    void g0p_sourceDensity(const PylithInt dim,
                           const PylithInt numS,
                           const PylithInt numA,
                           const PylithInt sOff[],
                           const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity

        // 3 + n
        const PylithInt i_sourceDensity = 3;

        const PylithInt _numS = 1; // Number passed on to g0p_source.

        const PylithInt numASource = 1; // Number passed on to g0p_source.
        const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
        const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

        pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                      NULL, NULL, NULL, NULL, NULL,
                                                      aOffSource, aOffSource_x, a, a_t, a_x,
                                                      t, x, numConstants, constants, g0);
    } // g0p_sourceDensity

    // ------------------------------------------------------------------------------
    /* g0p function for isotropic linear Poroelasticity plane strain with source density, gravity
     *
     */
    static inline
    void g0p_sourceDensity_grav(const PylithInt dim,
                                const PylithInt numS,
                                const PylithInt numA,
                                const PylithInt sOff[],
                                const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity

        // 2 + n
        const PylithInt i_sourceDensity = 4;

        const PylithInt _numS = 1; // Number passed on to g0p_source.

        const PylithInt numASource = 1; // Number passed on to g0p_source.
        const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
        const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

        pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                      NULL, NULL, NULL, NULL, NULL,
                                                      aOffSource, aOffSource_x, a, a_t, a_x,
                                                      t, x, numConstants, constants, g0);
    } // g0p_sourceDensity_grav

    // ------------------------------------------------------------------------------
    /*
     * g0p function for isotropic linear Poroelasticity plane strain with source density, and body force.
     */
    static inline
    void g0p_sourceDensity_body(const PylithInt dim,
                                const PylithInt numS,
                                const PylithInt numA,
                                const PylithInt sOff[],
                                const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity

        // 3 + n
        const PylithInt i_sourceDensity = 4;

        const PylithInt _numS = 1; // Number passed on to g0p_source.

        const PylithInt numASource = 1; // Number passed on to g0p_source.
        const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
        const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

        pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                      NULL, NULL, NULL, NULL, NULL,
                                                      aOffSource, aOffSource_x, a, a_t, a_x,
                                                      t, x, numConstants, constants, g0);
    } // g0p_sourceDensity_body

    // ------------------------------------------------------------------------------
    /*
     * g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
     */
    static inline
    void g0p_sourceDensity_grav_body(const PylithInt dim,
                                     const PylithInt numS,
                                     const PylithInt numA,
                                     const PylithInt sOff[],
                                     const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.

        // Poroelasticity

        // 3 + n
        const PylithInt i_sourceDensity = 5;

        const PylithInt _numS = 1; // Number passed on to g0p_source.

        const PylithInt numASource = 1; // Number passed on to g0p_source.
        const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
        const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

        pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                      NULL, NULL, NULL, NULL, NULL,
                                                      aOffSource, aOffSource_x, a, a_t, a_x,
                                                      t, x, numConstants, constants, g0);
    } // g0p_sourceDensity_grav_body

    /* -------------------------------------------------------------------------- */
    /*                           LHS Jacobian                                     */
    /* -------------------------------------------------------------------------- */

    // -----------------------------------------------------------------------------
    /*
     * Jg0ee - Jf0 function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        assert(aOff);
        assert(a);

        Jf0[0] = -1.0;
    } // Jg0ee

    // -----------------------------------------------------------------------------
    /*
     * Jf1eu - Jf1 function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf1[]) {
        for (PylithInt d = 0; d < dim; ++d) {
            Jf1[d*dim+d] = 1.0;
        } // for
    } // Jf1eu

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0vu function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0vu_implicit(const PylithInt dim,
                        const PylithInt numS,
                        const PylithInt numA,
                        const PylithInt sOff[],
                        const PylithInt sOff_x[],
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
        assert(aOff);
        assert(a);

        // Incoming auxiliary fields.

        for (PylithInt d = 0; d < dim; ++d) {
            Jf0[d * dim + d] += s_tshift;
        } // for
    } // Jf0vu_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0vv function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0vv_implicit(const PylithInt dim,
                        const PylithInt numS,
                        const PylithInt numA,
                        const PylithInt sOff[],
                        const PylithInt sOff_x[],
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
        assert(aOff);
        assert(a);

        // Incoming auxiliary fields.

        for (PylithInt d = 0; d < dim; ++d) {
            Jf0[d * dim + d] -= 1.0;
        } // for
    } // Jf0vv_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0vv function for poroelasticity equation, dynamic
     */
    static inline
    void Jf0vv_explicit(const PylithInt dim,
                        const PylithInt numS,
                        const PylithInt numA,
                        const PylithInt sOff[],
                        const PylithInt sOff_x[],
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
        // Incoming auxiliary fields.
        const PylithInt i_solid_density = 0;
        const PylithInt i_fluid_density = 1;
        const PylithInt i_porosity = 3;

        const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];

        for (PetscInt i = 0; i < dim; ++i) {
            Jf0[i*dim+i] += s_tshift * bulkDensity;
        } // for
    } // Jf0vv_explicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0pdotp function for poroelasticity equation, quasistatic.
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
        assert(aOff);
        assert(a);

        Jf0[0] += s_tshift;
    } // Jf0pdotp

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0pdotpdot function for poroelasticity equation, quasistatic.
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
        assert(aOff);
        assert(a);

        Jf0[0] -= 1.0;
    } // Jg0pdotpdot

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0edote function for poroelasticity equation, quasistatic.
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
        assert(aOff);
        assert(a);

        Jf0[0] += s_tshift;
    } // Jf0edote

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0edotedot function for poroelasticity equation, quasistatic.
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
        assert(aOff);
        assert(a);

        Jf0[0] -= 1.0;
    } // Jg0edotedot

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric strain or stress.
     *
     * Order of output components.
     *   2D: xx, yy, zz, xy
     *   3D: xx, yy, zz, xy, yz, xz
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void deviatoric(const pylith::fekernels::Tensor& tensor,
                    pylith::fekernels::Tensor* deviatoricTensor) {
        assert(deviatoricTensor);

        const PylithReal mean = (tensor.xx + tensor.yy + tensor.zz) / 3.0;
        deviatoricTensor->xx = tensor.xx - mean;
        deviatoricTensor->yy = tensor.yy - mean;
        deviatoricTensor->zz = tensor.zz - mean;
        deviatoricTensor->xy = tensor.xy;
        deviatoricTensor->yz = tensor.yz;
        deviatoricTensor->xz = tensor.xz;
    } // deviatoricStain

    // --------------------------------------------------------------------------------------------
    /** Calculate strain as a vector.
     *
     * Order of output components.
     *   2D: xx, yy, zz, xy
     *   3D: xx, yy, zz, xy, yz, xz
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void strain_asVector(const PylithInt dim,
                         const PylithInt numS,
                         const PylithInt sOff[],
                         const PylithInt sOff_x[],
                         const PylithScalar s[],
                         const PylithScalar s_t[],
                         const PylithScalar s_x[],
                         const PylithScalar x[],
                         strainfn_type strainFn,
                         const TensorOps& tensorOps,
                         PylithScalar strainVector[]) {
        assert(strainVector);

        Tensor strain;
        strainFn(dim, numS, sOff, sOff_x, s, s_t, s_x, x, &strain);
        tensorOps.toVector(strain, strainVector);
    } // infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Calculate stress as a vector.
     *
     * Order of output components.
     *   2D: xx, yy, zz, xy
     *   3D: xx, yy, zz, xy, yz, xz
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void stress_asVector(const PylithInt dim,
                         const PylithInt numS,
                         const PylithInt numA,
                         const PylithInt sOff[],
                         const PylithInt sOff_x[],
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
                         strainfn_type strainFn,
                         stressfn_type stressFn,
                         const TensorOps& tensorOps,
                         PylithScalar stressVector[]) {
        assert(stressVector);

        Tensor strain;
        strainFn(dim, numS, sOff, sOff_x, s, s_t, s_x, x, &strain);

        Tensor stress;
        stressFn(dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
                 t, x, numConstants, constants, strain, tensorOps, &stress);

        tensorOps.toVector(stress, stressVector);
    } // cauchyStress_asVector

}; // Poroelasticity

// ---------------------------------------------------------------------------------------------------------------------

/// Kernels specific to poroelasticity plane strain.
class pylith::fekernels::PoroelasticityPlaneStrain {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Calculate infinitesimal strain tensor for 2-D plane strain poroelasticity.
     *
     * * Order of output components is xx, xy, yx, yy.
     * Order of output components is xx, yy, zz, xy.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithScalar x[],
                             pylith::fekernels::Tensor* strain) {
        const PylithInt _dim = 2;

        assert(_dim == dim);
        assert(numS >= 1);
        assert(sOff_x);
        assert(s_x);
        assert(strain);

        // Incoming solution field.
        const PylithInt i_disp = 0;
        const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

        strain->xx = disp_x[0*_dim+0];
        strain->yy = disp_x[1*_dim+1];
        strain->zz = 0.0;
        strain->xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);
        strain->yz = 0.0;
        strain->xz = 0.0;
    } // infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with infinitesimal strain for 2D plane strain poroelasticity.
     *
     * Order of output components is xx, yy, zz, xy.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain_asVector(const PylithInt dim,
                                      const PylithInt numS,
                                      const PylithInt numA,
                                      const PylithInt sOff[],
                                      const PylithInt sOff_x[],
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
                                      PylithScalar strainVector[]) {
        const PylithInt _dim = 2;
        assert(_dim == dim);

        Poroelasticity::strain_asVector(_dim, numS, sOff, sOff_x, s, s_t, s_x, x,
                                        infinitesimalStrain, Tensor::ops2D, strainVector);
    } // infinitesimalStrain_asVector3D

}; // PoroelasticityPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to poroelasticity in 3D.
class pylith::fekernels::Poroelasticity3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Calculate infinitesimal strain for 3D poroelasticity.
     *
     * Order of output components is xx, xy, xz, yx, yy, yz, zx, zy, zz.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain(const PylithInt dim,
                             const PylithInt numS,
                             const PylithInt sOff[],
                             const PylithInt sOff_x[],
                             const PylithScalar s[],
                             const PylithScalar s_t[],
                             const PylithScalar s_x[],
                             const PylithScalar x[],
                             pylith::fekernels::Tensor* strain) {
        const PylithInt _dim = 3;

        assert(_dim == dim);
        assert(numS >= 1);
        assert(sOff_x);
        assert(s_x);
        assert(strain);

        // Incoming solution field.
        const PylithInt i_disp = 0;
        const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

        strain->xx = disp_x[0*_dim+0];
        strain->yy = disp_x[1*_dim+1];
        strain->zz = disp_x[2*_dim+2];
        strain->xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);
        strain->yz = 0.5*(disp_x[1*_dim+2] + disp_x[2*_dim+1]);
        strain->xz = 0.5*(disp_x[0*_dim+2] + disp_x[2*_dim+0]);
    } // infinitesimalStrain

    // --------------------------------------------------------------------------------------------
    /** Calculate vector with infinitesimal strain for 3D poroelasticity.
     *
     * Order of output components is xx, yy, zz, xy, yz, xz.
     *
     * Solution fields: [disp(dim)]
     */
    static inline
    void infinitesimalStrain_asVector(const PylithInt dim,
                                      const PylithInt numS,
                                      const PylithInt numA,
                                      const PylithInt sOff[],
                                      const PylithInt sOff_x[],
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
                                      PylithScalar strainVector[]) {
        const PylithInt _dim = 3;
        assert(_dim == dim);

        Poroelasticity::strain_asVector(_dim, numS, sOff, sOff_x, s, s_t, s_x, x,
                                        infinitesimalStrain, Tensor::ops3D, strainVector);
    } // infinitesimalStrain_asVector

}; // Poroelasticity3D

#endif // pylith_fekernels_poroelasticity_hh

// End of file
