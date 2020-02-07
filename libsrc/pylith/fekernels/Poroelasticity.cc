/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/Poroelasticity.hh"

#include <cassert> // USES assert()
#include <iostream> // use to output data to screen

// =====================================================================================================================
// Generic poroelasticity kernels for inertia and body forces.
// =====================================================================================================================

/* -------------------------------------------------------------------------- */
/*                           LHS Residuals                                    */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/*                           RHS Residuals                                    */
/* -------------------------------------------------------------------------- */
// Quasi-Static

// =============================================================================
// Displacement
// =============================================================================
// ---------------------------------------------------------------------------------------------------------------------
// g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
void
pylith::fekernels::Poroelasticity::g0v_grav(const PylithInt dim,
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
    const PylithInt i_porosity     = 0;
    const PylithInt i_density      = 1;
    const PylithInt i_fluidDensity = 2;

    const PylithInt i_gravityField = 4;

    // assert(_numS == numS);
    // assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithScalar density = (1 - a[aOff[i_porosity]]) * a[aOff[i_density]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density * gravityField[i];
    } // for

} // g0v_grav

// ---------------------------------------------------------------------------------------------------------------------
// g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
void
pylith::fekernels::Poroelasticity::g0v_bodyforce(const PylithInt dim,
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
  const PylithInt i_bodyForce = 4;
  assert(aOff);
  assert(aOff[i_bodyForce] >= 0);
  assert(a);

  const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

  for (PylithInt i = 0; i < dim; ++i) {
    g0[i] += bodyForce[i];
  } // for
} // g0v_bodyforce


// ----------------------------------------------------------------------
//g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
void
pylith::fekernels::Poroelasticity::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt i_porosity = 0;
    const PylithInt i_density = 1;
    const PylithInt i_fluidDensity = 2;
    const PylithInt i_gravityField = 4;
    const PylithInt i_bodyForce = 5;

    assert(aOff);
    assert(aOff_x);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithScalar density = (1 - a[aOff[i_porosity]]) * a[aOff[i_density]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    // gravity field
    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density * gravityField[i];
    } // for

    // body force
    for (PylithInt i = 0; i < dim; ++i) {
      g0[i] += bodyForce[i];
    } // for

} // g0v_gravbodyforce

// =============================================================================
// Pressure
// =============================================================================

// ----------------------------------------------------------------------
//g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity(const PylithInt dim,
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
                                             PylithScalar g0p[]) {
    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 4;

    assert(aOff);
    assert(aOff[i_sourceDensity] >= 0);
    assert(a);

    const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0p[i] += sourceDensity[i];
    } // for
} // g0p_source

// ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body(const PylithInt dim,
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
                                                                       PylithScalar g0p[]) {


    //const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 4;

    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_sourceDensity(dim, _numS, numASource,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffSource, aOffSource_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0p);
} // g0p_sourceDensity_grav_body



// =============================================================================
// Volumetric Strain
// =============================================================================
// ----------------------------------------------------------------------
// g0E function for isotropic linear Poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::g0e_trace_strain(const PylithInt dim,
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
                                                   PylithScalar g0E[]) {
    //const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_trace = 2;

    // Incoming auxiliary fields.

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffTrace[2] = { sOff[i_disp], sOff[i_trace] };
    const PylithInt sOffTrace_x[2] = { sOff_x[i_disp], sOff_x[i_trace] };

    pylith::fekernels::Poroelasticity::trace_strainCal(dim, _numS, 0,
                                                         sOffTrace, sOffTrace_x, s, s_t, s_x,
                                                         NULL, NULL, NULL, NULL, NULL,
                                                         t, x, numConstants, constants, g0E);
} // g0e_trace_strain

// ----------------------------------------------------------------------
/*
 *
 */
void
pylith::fekernels::Poroelasticity::trace_strainCal(const PylithInt dim,
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
                                                     PylithScalar g0E[]) {


    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_trace = 1;

    const PylithScalar* disp = &s[sOff[i_disp]];
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar trace = s[sOff[i_trace]];

    PylithScalar strainTrace = 0.0;

    if (dim == 1) {
        strainTrace = disp_x[0*dim+0];
    } else if (dim == 2) {
        strainTrace = disp_x[0*dim+0] + disp_x[1*dim+1];
    } else if (dim == 3) {
        strainTrace = disp_x[0*dim+0] + disp_x[1*dim+1] + disp_x[2*dim+2];
    } //elseif

    // for (PylithInt i = 0; i < dim; ++i) {
    //    g0E[i] += strainTrace - trace;
    //}

    g0E[0] += strainTrace - trace;


} // trace_strainCal



/* -------------------------------------------------------------------------- */
/*                           LHS Jacobian                                     */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/*                           RHS Jacobian                                     */
/* -------------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
//Jg0ee - Jg0 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jg0ee(const PylithInt dim,
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
                                                PylithScalar Jg0[]) {

    assert(aOff);
    assert(a);

    Jg0[0] += -1;
} // Jg0ee

// -----------------------------------------------------------------------------
// Jg1eu - Jg1 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jg1eu(const PylithInt dim,
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
                                                PylithScalar Jg1[]) {
    //const PylithInt _dim = 2;
    PylithInt i;
    assert(aOff);
    assert(a);

    for (i = 0; i < dim; ++i) {
        Jg1[i*dim+i] += 1.0;
    } // for
} // Jg1eu

// =====================================================================================================================
// Kernels for linear isotropic poroelasticity
// =====================================================================================================================

// ---------------------------------------------------------------------------------------------------------------------
/* Calculate Cauchy strain for 2-D plane strain elasticity.
 *
 * Order of output components is xx, yy, zz, xy for 2D,
 * Order of output components is xx, yy, zz, xy, yz, xz for 3D.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::Poroelasticity::cauchyStrain(const PylithInt dim,
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
                                                       PylithScalar strain[]) {
//    const PylithInt _dim = 2;

//    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    if (dim == 2) {
        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_zz = 0.0;
        const PylithScalar strain_xy = 0.5*(disp_x[0*dim+1] + disp_x[1*dim+0]);

        strain[0] = strain_xx;
        strain[1] = strain_yy;
        strain[2] = strain_zz;
        strain[3] = strain_xy;

    } else if (dim == 3) {

        const PylithScalar strain_xx = disp_x[0*dim+0];
        const PylithScalar strain_yy = disp_x[1*dim+1];
        const PylithScalar strain_zz = disp_x[2*dim+2];
        const PylithScalar strain_xy = 0.5*(disp_x[0*dim+1] + disp_x[1*dim+0]);
        const PylithScalar strain_yz = 0.5*(disp_x[1*dim+2] + disp_x[2*dim+1]);
        const PylithScalar strain_xz = 0.5*(disp_x[0*dim+2] + disp_x[2*dim+0]);

        strain[0] = strain_xx;
        strain[1] = strain_yy;
        strain[2] = strain_zz;
        strain[3] = strain_xy;
        strain[4] = strain_yz;
        strain[5] = strain_xz;
    }
} // cauchyStrain
