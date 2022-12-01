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

#include "pylith/fekernels/FaultCohesiveKinPoro.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for prescribed fault slip.
 *
 * Solution fields: [disp(dim), vel(dim, if dynamic), pressure(1), trace_strain(1), lagrange(dim), fault_pressure(1)]
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 *
 ***** Required fields
 * thickness(1),                       0
 * porosity(1),                        1
 * beta_p(1),                          2
 * beta_sigma(1),                      3
 * permeability_tangential(1),         4
 * permeability_normal(1),             5
 * fluid_viscosity(1),                 6+
 * body_force(dim),                    numA - 3
 * source (1),                         numA - 2
 * slip(dim)                           numA - 1
 * ** TO DO **
 *
 * Solution fields
 *  ***** Required fields
 * slip(3),                            0        (soff not implemented, trivial)
 * velocity(3),                        1
 * pressure(1),                        numS - 4 (quasi-static) numS - 3 (dynamic)
 * trace_strain [if quasi-static](1),  numS - 3
 * lagrange,                           numS - 2
 * fault_pressure,                     numS - 1
 * ** TO DO **
 * (**NEED IMPLEMENTATION**) *
 * (**NEED IMPLEMENTATION**) *
 * ======================================================================
 */

namespace pylith {
    namespace fekernels {
        class _FaultCohesiveKinPoro {
public:

            /** Get offset in s where velocity subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of velocity subfield in s.
             */
            static PylithInt velocity_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where pressure subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of pressure subfield in s.
             */
            static PylithInt pressure_sOff_quasistatic(const PylithInt sOff[],
                                                       const PylithInt numS);

            /** Get offset in s where pressure subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of pressure subfield in s.
             */
            static PylithInt pressure_sOff_dynamic(const PylithInt sOff[],
                                                   const PylithInt numS);

            /** Get offset in s where trace_strain subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the velocity
             * subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of trace_strain subfield in s.
             */
            static PylithInt trace_strain_sOff(const PylithInt sOff[],
                                               const PylithInt numS);

            /** Get offset in s where Lagrange multiplier subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of Lagrange multiplier subfield in s.
             */
            static PylithInt lagrange_sOff(const PylithInt sOff[],
                                           const PylithInt numS);

            /** Get offset in s where fault_pressure subfield starts.
             *
             * Normally this would be sOff, but sOff doesn't account for having DOF for the two sides of the fault
             * passed to the hybrid kernels. This functions computes the correct offset into s for the Lagrange
             * multiplier subfield.
             *
             * @param[in] sOff Offset of registered subfields in solution field [numS].
             * @param[in] numS Number of registered subfields in solution field.
             *
             * @returns Offset of fault_pressure subfield in s.
             */
            static PylithInt fault_pressure_sOff(const PylithInt sOff[],
                                                 const PylithInt numS);

            /* Compute tangential directions for 3-D fault.
             *
             * @param[in] dim Spatial dimension.
             * @param[in] refDir1 First choice for reference direction.
             * @param[in] refDir2 Second choice for reference direction if first fails.
             * @param[in] normDir Normal direction.
             * @param[out] tanDir1 First tangential direction.
             * @param[out] tanDIr2 Second tangential direction.
             */
            static void tangential_directions(const PylithInt dim,
                                              const PylithScalar refDir1[],
                                              const PylithScalar refDir2[],
                                              const PylithScalar normDir[],
                                              PylithScalar tanDir1[],
                                              PylithScalar tanDir2[]);

        }; // _FaultCohesiveKinPoro
    } // fekernels
} // pylith

// ----------------------------------------------------------------------
// Get offset in s where velocity subfield starts.
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::velocity_sOff(const PylithInt sOff[],
                                                        const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = 1; // [displacement, velocity, ...]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // velocity_sOff


// Get offset in s where pressure subfield starts.
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::pressure_sOff_quasistatic(const PylithInt sOff[],
                                                                    const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 4; // [displacement, velocity, pressure...]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // pressure_sOff_quasistatic


// Get offset in s where pressure subfield starts.
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::pressure_sOff_dynamic(const PylithInt sOff[],
                                                                const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 3; // [displacement, velocity, pressure...]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // pressure_sOff_dynamic


// Get offset in s where trace_strain subfield starts.
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::trace_strain_sOff(const PylithInt sOff[],
                                                            const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 3; // [displacement, velocity, pressure, trace_strain...]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // trace_strain_sOff


// ----------------------------------------------------------------------
// Get offset in s where Lagrange multiplier field starts.
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::lagrange_sOff(const PylithInt sOff[],
                                                        const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 2; // Don't include last 2 field (Lagrange multiplier, fault_pressure)
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for
    return off;
} // lagrange_sOff


// ----------------------------------------------------------------------
// Get offset in s where pressure_fault subfield starts.
// Seems that the offset contributed by Lagrange multiplier does not need *2
// ** TO DO **
// ** NEEDS VERFICATION **
PylithInt
pylith::fekernels::_FaultCohesiveKinPoro::fault_pressure_sOff(const PylithInt sOff[],
                                                              const PylithInt numS) {
    PylithInt off = 0;
    const PylithInt numCount = numS - 2; // [..., Lagrange multiplier, fault_pressure]
    for (PylithInt i = 0; i < numCount; ++i) {
        off += 2 * (sOff[i + 1] - sOff[i]);
    } // for

    off += (sOff[numS - 1] - sOff[numS - 2]);
    return off;
} // fault_pressure_sOff


// ----------------------------------------------------------------------
// Compute tangential directions from reference direction (first and second choice) and normal direction in 3-D.
void
pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(const PylithInt dim,
                                                                const PylithScalar refDir1[],
                                                                const PylithScalar refDir2[],
                                                                const PylithScalar normDir[],
                                                                PylithScalar tanDir1[],
                                                                PylithScalar tanDir2[]) {
    assert(3 == dim);
    assert(refDir1);
    assert(refDir2);
    assert(normDir);
    assert(tanDir1);
    assert(tanDir2);

    const PylithInt _dim = 3;
    PylithScalar refDir[3] = {refDir1[0], refDir1[1], refDir1[2]};
    if (fabs(refDir[0] * normDir[0] + refDir[1] * normDir[1] + refDir[2] * normDir[2]) > 0.98) {
        for (PylithInt i = 0; i < _dim; ++i) {
            refDir[i] = refDir2[i];
        } // for
    } // if

    // refDir x normDir, normalization required
    tanDir1[0] = +refDir[1] * normDir[2] - refDir[2] * normDir[1];
    tanDir1[1] = +refDir[2] * normDir[0] - refDir[0] * normDir[2];
    tanDir1[2] = +refDir[0] * normDir[1] - refDir[1] * normDir[0];

    PylithScalar _norm = sqrt(tanDir1[0] * tanDir1[0] + tanDir1[1] * tanDir1[1] + tanDir1[2] * tanDir1[2]);
    tanDir1[0] /= _norm;
    tanDir1[1] /= _norm;
    tanDir1[2] /= _norm;

    // normDir x tanDir1, normalization not required
    tanDir2[0] = +normDir[1] * tanDir1[2] - normDir[2] * tanDir1[1];
    tanDir2[1] = +normDir[2] * tanDir1[0] - normDir[0] * tanDir1[2];
    tanDir2[2] = +normDir[0] * tanDir1[1] - normDir[1] * tanDir1[0];
} // _tangential_directions


// ----------------------------------------------------------------------
// Residual Functions
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// f0 function for elasticity equation: f0u = -\lambda (pos side), +\lambda (neg side).
void
pylith::fekernels::FaultCohesiveKinPoro::f0u_neg(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt fOffN = 0;
    const PylithInt sOffLagrange = sOff[numS - 2];
    const PylithScalar *lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffN + i] += lagrange[i];
        if (f0[fOffN + i] != f0[fOffN + i]) {
            PetscPrintf(PETSC_COMM_WORLD, "Error in f0u_neg \n");
        }
    } // for
} // f0u_neg


// ----------------------------------------------------------------------
// f0 function for poroelasticity equation: f0u = -\lambda (pos side), +\lambda (neg side).
void
pylith::fekernels::FaultCohesiveKinPoro::f0u_pos(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields

    const PylithInt fOffP = 0;
    const PylithInt sOffLagrange = sOff[numS - 2];
    const PylithScalar *lagrange = &s[sOffLagrange];

    for (PylithInt i = 0; i < spaceDim; ++i) {
        f0[fOffP + i] += -lagrange[i];
        if (f0[fOffP + i] != f0[fOffP + i]) {
            PetscPrintf(PETSC_COMM_WORLD, "Error in f0u_pos \n");
        }
    } // for
} // f0u_pos


// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_neg(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = numS - 1;

    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt sOffPressureN = sOff[i_pressure];
    const PylithInt sOffPressureFault = sOff[i_fault_pressure];
    const PylithInt fOffN = 0;

    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar pressureN = s[sOffPressureN];
    // const PylithScalar pressureP = s[sOffPressureP];
    const PylithScalar pressureFault = s[sOffPressureFault];

    f0[fOffN] += 0.0;
    if (f0[fOffN] != f0[fOffN]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_neg \n");
    }
    // f0[fOffN] += permeabilityNormal / fluidViscosity *
    //              ((pressureN - pressureFault) / thickness);
    // f0[fOffP] += permeabilityNormal / fluidViscosity *
    //             ((pressureP - pressureFault) / thickness);
} // f0p_neg


// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_pos(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = 4;

    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt sOffPressureP = sOff[i_pressure];
    const PylithInt sOffPressureFault = sOff[i_fault_pressure];

    const PylithInt fOffP = 0;

    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    // const PylithScalar pressureN = s[sOffPressureN];
    const PylithScalar pressureP = s[sOffPressureP];
    const PylithScalar pressureFault = s[sOffPressureFault];

    f0[fOffP] += 0.0;
    if (f0[fOffP] != f0[fOffP]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_pos \n");
    }
    // f0[fOffN] += permeabilityNormal / fluidViscosity *
    //             ((pressureN - pressureFault) / thickness);
    // f0[fOffP] += permeabilityNormal / fluidViscosity *
    //              ((pressureP - pressureFault) / thickness);
} // f0p_pos


// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_body_neg(const PylithInt dim,
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
                                                      const PylithReal n[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = 4;

    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;

    const PylithInt fOffN = 0;

    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_body_force]];
    const PylithScalar pressureN = s[sOff[i_pressure]];
    const PylithScalar pressureFault = s[sOff[i_fault_pressure]];

    PylithScalar nDotBodyForce = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotBodyForce += n[i] * bodyForce[i];
    }

    f0[fOffN] += 0.0;

    if (f0[fOffN] != f0[fOffN]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_body_neg \n");
    }

    // f0[fOffN] += permeabilityNormal / fluidViscosity *
    //              ((pressureN - pressureFault) / thickness + nDotBodyForce);
    // f0[fOffP] += permeabilityNormal / fluidViscosity *
    //             ((pressureP - pressureFault) / thickness - nDotBodyForce);
} // f0p_body_neg


// ----------------------------------------------------------------------
// f0 function for bulk pressure: f0p = [\kappa_{cz} / \mu * ((p^+ - p^f)/h - n \cdot f_f),
//                                       \kappa_{cz} / \mu * ((p^- - p^f)/h + n \cdot f_f)]
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_body_pos(const PylithInt dim,
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
                                                      const PylithReal n[],
                                                      const PylithInt numConstants,
                                                      const PylithScalar constants[],
                                                      PylithScalar f0[]) {
    assert(sOff);
    assert(s);
    assert(f0);

    assert(numS >= 5);
    assert(a);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = 4;

    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;

    const PylithInt fOffP = 0;

    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce = &a[aOff[i_body_force]];
    const PylithScalar pressureP = s[sOff[i_pressure]];
    const PylithScalar pressureFault = s[sOff[i_fault_pressure]];

    PylithScalar nDotBodyForce = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotBodyForce += n[i] * bodyForce[i];
    }

    f0[fOffP] += 0.0;

    if (f0[fOffP] != f0[fOffP]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_body_pos \n");
    }

    // f0[fOffN] += permeabilityNormal / fluidViscosity *
    //             ((pressureN - pressureFault) / thickness + nDotBodyForce);
    // f0[fOffP] += permeabilityNormal / fluidViscosity *
    //              ((pressureP - pressureFault) / thickness - nDotBodyForce);
} // f0p_body_pos


// ----------------------------------------------------------------------
// f0 function for slip constraint equation: f0\lambda = (u^+ - u^-) - d
void
pylith::fekernels::FaultCohesiveKinPoro::f0l_u(const PylithInt dim,
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
                                               const PylithReal n[],
                                               const PylithInt numConstants,
                                               const PylithScalar constants[],
                                               PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slip = numA - 1;
    const PylithInt i_disp = 0;

    const PylithScalar *slip = &a[aOff[i_slip]];

    const PylithInt sOffDispN = sOff[i_disp];
    const PylithInt sOffDispP = sOffDispN + spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar *dispN = &s[sOffDispN];
    const PylithScalar *dispP = &s[sOffDispP];

    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {n[1], -n[0]};
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXY = n[i] * slip[0] + tanDir[i] * slip[1];
            f0[fOffLagrange + i] += -dispP[i] + dispN[i] + slipXY;
            if (f0[fOffLagrange + i] != f0[fOffLagrange + i]) {
                PetscPrintf(PETSC_COMM_WORLD, "Error in f0l_u \n");
            }
        } // for
        break;
    } // case 2
    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipXYZ = n[i] * slip[0] + tanDir1[i] * slip[1] + tanDir2[i] * slip[2];
            f0[fOffLagrange + i] += -dispP[i] + dispN[i] + slipXYZ;
            if (f0[fOffLagrange + i] != f0[fOffLagrange + i]) {
                PetscPrintf(PETSC_COMM_WORLD, "Error in f0l_u \n");
            }
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_u


// ----------------------------------------------------------------------
// f0 function for slip rate constraint equation: f0\lambda = (v^+ - v^-) - \dot{d}
void
pylith::fekernels::FaultCohesiveKinPoro::f0l_v(const PylithInt dim,
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
                                               const PylithReal n[],
                                               const PylithInt numConstants,
                                               const PylithScalar constants[],
                                               PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(s);
    assert(a);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 6);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_slipRate = numA - 1;

    const PylithScalar *slipRate = &a[aOff[i_slipRate]];

    const PylithInt sOffVelN = _FaultCohesiveKinPoro::velocity_sOff(sOff, numS);
    const PylithInt sOffVelP = sOffVelN + spaceDim;
    const PylithInt fOffLagrange = 0;

    const PylithScalar *velN = &s[sOffVelN];
    const PylithScalar *velP = &s[sOffVelP];

    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipRateXY = n[i] * slipRate[0] + tanDir[i] * slipRate[1];
            f0[fOffLagrange + i] += -velP[i] + velN[i] + slipRateXY;
            if (f0[fOffLagrange + i] != f0[fOffLagrange + i]) {
                PetscPrintf(PETSC_COMM_WORLD, "Error in f0l_v \n");
            }
        } // for
        break;
    } // case 2
    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);

        for (PylithInt i = 0; i < _spaceDim; ++i) {
            const PylithScalar slipRateXYZ = n[i] * slipRate[0] + tanDir1[i] * slipRate[1] + tanDir2[i] * slipRate[2];
            f0[fOffLagrange + i] += -velP[i] + velN[i] + slipRateXYZ;
            if (f0[fOffLagrange + i] != f0[fOffLagrange + i]) {
                PetscPrintf(PETSC_COMM_WORLD, "Error in f0l_v \n");
            }
        } // for
        break;
    } // case 3
    default:
        assert(0);
    } // switch
} // f0l_v


// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (-n \cdot \lambda ))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_fault(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_lagrange = 3;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_fault_permeabilility = 4;
    const PylithInt i_fluid_viscosity = 5;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeabilility]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    // // Hardcoded test values
    // const PylithScalar porosity = 0.1;
    // const PylithScalar betaP = 1.0;
    // const PylithScalar betaSigma = 1.0;
    // const PylithScalar faultPermeability[4] = {1.0, 1.0, 0.0, 0.0};

    // Pressure and pressure_t
    const PylithInt sOffpressureN = sOff[i_pressure];
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffLagrange = sOff[i_lagrange];
    const PylithInt sOffpressureFault = sOff[i_fault_pressure];

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    // const PylithScalar *faultLagrange_t = &s_t[sOffLagrange];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];
    /**
     *  // Trace_strain, no transformation required
     *
     *  const PylithInt sOffTraceStrainN = _FaultCohesiveKinPoro::trace_strain_sOff(sOff, numS);
     *  const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
     *  const PylithScalar traceStrainN = s[sOffTraceStrainN];
     *  const PylithScalar traceStrainP = s[sOffTraceStrainP];
     *
     *  // \sigma_nn, requires transformation
     *  // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
     *  // strain_nn = n_i u_{i,j} n_j
     *  PylithScalar strain_nnN = 0.;
     *  PylithScalar strain_nnP = 0.;
     *  for (PylithInt i = 0; i < spaceDim; ++i)
     *  {
     *      for (PylithInt j = 0; j < spaceDim; ++j)
     *      {
     *          strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
     *          strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
     *      }
     *  }
     *
     *  const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN *
     * strain_nnN;
     *  const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP *
     * strain_nnP;
     */
    const PylithScalar *faultLagrange_t = &s_t[sOffLagrange];
    PylithScalar nDotLagrange_t = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotLagrange_t -= n[i] * faultLagrange_t[i];
    }

    // f0[fOffp_fault] += porosity * (betaP * (pressureP_t + 2. * pressureFault_t + pressureN_t) / 4. + betaSigma *
    // nDotLagrange_t);
    f0[fOffp_fault] += porosity * (betaP * pressureFault_t + betaSigma * nDotLagrange_t);
    if (f0[fOffp_fault] != f0[fOffp_fault]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_fault \n");
    }

} // f0p_fault


// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_fault_body(const PylithInt dim,
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
                                                        const PylithReal n[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_lagrange = 3;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_fault_permeabilility = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeabilility]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];

    // Divergence of body force, should be naturally zero in direction n,
    // but here still correct it;
    PylithScalar bodyForce_div = 0.;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        bodyForce_div += bodyForce_x[i * spaceDim + i];
    }

    for (PylithInt i = 0; i < spaceDim; ++i) {
        for (PylithInt j = 0; j < spaceDim; ++j) {
            bodyForce_div -= n[i] * n[j] * bodyForce_x[i * spaceDim + j];
        }
    }

    // Pressure and pressure_t
    const PylithInt sOffpressureN = sOff[i_pressure];
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffLagrange = sOff[i_lagrange];
    const PylithInt sOffpressureFault = sOff[i_fault_pressure];

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];

    // Trace_strain, no transformation required
    /**
     *  const PylithInt sOffTraceStrainN = _FaultCohesiveKinPoro::trace_strain_sOff(sOff, numS);
     *  const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
     *  const PylithScalar traceStrainN = s[sOffTraceStrainN];
     *  const PylithScalar traceStrainP = s[sOffTraceStrainP];
     *
     *  // \sigma_nn, requires transformation
     *  // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
     *  // strain_nn = n_i u_{i,j} n_j
     *  PylithScalar strain_nnN = 0.;
     *  PylithScalar strain_nnP = 0.;
     *  for (PylithInt i = 0; i < spaceDim; ++i)
     *  {
     *      for (PylithInt j = 0; j < spaceDim; ++j)
     *      {
     *          strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
     *          strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
     *      }
     *  }
     *
     *  const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN *
     * strain_nnN;
     *  const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP *
     * strain_nnP;
     */
    const PylithScalar *faultLagrange_t = &s_t[sOffLagrange];
    PylithScalar nDotLagrange_t = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotLagrange_t -= n[i] * faultLagrange_t[i];
    }

    f0[fOffp_fault] += porosity * (betaP * (pressureP_t + 2. * pressureFault_t + pressureN_t) / 4. +
                                   betaSigma * nDotLagrange_t);
    if (f0[fOffp_fault] != f0[fOffp_fault]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_fault_body \n");
    }

} // f0p_fault_body


// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_fault_source(const PylithInt dim,
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
                                                          const PylithReal n[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_lagrange = 3;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;
    const PylithInt i_source = numA - 2;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeability]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];
    const PylithScalar source = a[aOff[i_source]];

    // Pressure and pressure_t
    const PylithInt sOffpressureN = sOff[i_pressure];
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffLagrange = sOff[i_lagrange];
    const PylithInt sOffpressureFault = sOff[i_fault_pressure];

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    // const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    // const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    // const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    // const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];
    /**
     *  // Trace_strain, no transformation required
     *
     *  const PylithInt sOffTraceStrainN = _FaultCohesiveKinPoro::trace_strain_sOff(sOff, numS);
     *  const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
     *  const PylithScalar traceStrainN = s[sOffTraceStrainN];
     *  const PylithScalar traceStrainP = s[sOffTraceStrainP];
     *
     *  // \sigma_nn, requires transformation
     *  // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
     *  // strain_nn = n_i u_{i,j} n_j
     *  PylithScalar strain_nnN = 0.;
     *  PylithScalar strain_nnP = 0.;
     *  for (PylithInt i = 0; i < spaceDim; ++i)
     *  {
     *      for (PylithInt j = 0; j < spaceDim; ++j)
     *      {
     *          strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
     *          strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
     *      }
     *  }
     *
     *  const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN *
     * strain_nnN;
     *  const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP *
     * strain_nnP;
     */
    const PylithScalar *faultLagrange_t = &s_t[sOffLagrange];
    PylithScalar nDotLagrange_t = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotLagrange_t -= n[i] * faultLagrange_t[i];
    }

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * nDotLagrange_t) -
                       source;
    if (f0[fOffp_fault] != f0[fOffp_fault]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_fault_source \n");
    }

} // f0p_fault_source


// ----------------------------------------------------------------------
// f0 function for p_fault constraint equation:
// f0p_fault = porosity * (\beta^p * (\dot{p}^+ + 2 \dot{p}^f + \dot{p}^-)/4
//                         + \beta^\sigma * (\dot{\sigma}^+_{nn} + \dot{\sigma}^-_{nn}))
//             + \kappa_{fx} / \mu * \vnabla(2D) \cdot body_force
//             - \kappa_{fz} / \mu * (p^+ - 2p_f + p^-) / h^2 - source
void
pylith::fekernels::FaultCohesiveKinPoro::f0p_fault_body_source(const PylithInt dim,
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
                                                               const PylithReal n[],
                                                               const PylithInt numConstants,
                                                               const PylithScalar constants[],
                                                               PylithScalar f0[]) {
    assert(sOff);
    assert(aOff);
    assert(sOff_x);
    assert(a_x);
    assert(aOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(s_t);
    assert(f0);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
    const PylithInt i_disp_x = 0;

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_lagrange = 3;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;
    const PylithInt i_beta_sigma = 3;
    const PylithInt i_fault_permeabilility = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;
    const PylithInt i_source = numA - 2;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];
    const PylithScalar *faultPermeability = &a[aOff[i_fault_permeabilility]];
    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *bodyForce_x = &a_x[aOff_x[i_body_force]];
    const PylithScalar source = a[aOff[i_source]];

    // Divergence of body force, corrected by subtracting the nn term.
    PylithScalar bodyForce_div = 0.;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        bodyForce_div += bodyForce_x[i * spaceDim + i];
    }

    for (PylithInt i = 0; i < spaceDim; ++i) {
        for (PylithInt j = 0; j < spaceDim; ++j) {
            bodyForce_div -= n[i] * n[j] * bodyForce_x[i * spaceDim + j];
        }
    }

    // Pressure and pressure_t
    const PylithInt sOffpressureN = sOff[i_pressure];
    const PylithInt sOffpressureP = sOffpressureN + 1;
    const PylithInt sOffLagrange = sOff[i_lagrange];
    const PylithInt sOffpressureFault = sOff[i_fault_pressure];

    const PylithScalar pressureN = s[sOffpressureN];
    const PylithScalar pressureP = s[sOffpressureP];
    const PylithScalar pressureFault = s[sOffpressureFault];
    const PylithScalar pressureN_t = s_t[sOffpressureN];
    const PylithScalar pressureP_t = s_t[sOffpressureP];
    const PylithScalar pressureFault_t = s_t[sOffpressureFault];

    // ** TO DO **
    // Pull out drained bulk modulus and shear modulus from the surrounding bulks
    // const PylithScalar bulkModulusN = a[aOff[i_bulk_modulus_negative]];
    // const PylithScalar bulkModulusP = a[aOff[i_bulk_modulus_positive]];
    // const PylithScalar shearModulusN = a[aOff[i_shear_modulus_negative]];
    // const PylithScalar shearModulusP = a[aOff[i_shear_modulus_positive]];

    // Strain components
    const PylithInt sOffDispN_x = sOff_x[i_disp_x];
    const PylithInt sOffDispP_x = sOffDispN_x + spaceDim ^ 2;
    const PylithInt fOffp_fault = 0;

    const PylithScalar *dispN_x = &s_x[sOffDispN_x];
    const PylithScalar *dispP_x = &s_x[sOffDispP_x];
    /**
     *  // Trace_strain, no transformation required
     *
     *  const PylithInt sOffTraceStrainN = _FaultCohesiveKinPoro::trace_strain_sOff(sOff, numS);
     *  const PylithInt sOffTraceStrainP = sOffTraceStrainN + 1;
     *  const PylithScalar traceStrainN = s[sOffTraceStrainN];
     *  const PylithScalar traceStrainP = s[sOffTraceStrainP];
     *
     *  // \sigma_nn, requires transformation
     *  // \sigma_nn = (K_u - 2G/3) \epsilon_v + 2G strain_nn
     *  // strain_nn = n_i u_{i,j} n_j
     *  PylithScalar strain_nnN = 0.;
     *  PylithScalar strain_nnP = 0.;
     *  for (PylithInt i = 0; i < spaceDim; ++i)
     *  {
     *      for (PylithInt j = 0; j < spaceDim; ++j)
     *      {
     *          strain_nnN += n[i] * dispN_x[i * spaceDim + j] * n[j];
     *          strain_nnP += n[i] * dispP_x[i * spaceDim + j] * n[j];
     *      }
     *  }
     *
     *  const PylithScalar stress_nnN = (bulkModulusN - 2. * shearModulusN / 3.) * traceStrainN + 2. * shearModulusN *
     * strain_nnN;
     *  const PylithScalar stress_nnP = (bulkModulusP - 2. * shearModulusP / 3.) * traceStrainP + 2. * shearModulusP *
     * strain_nnP;
     */

    const PylithScalar *faultLagrange_t = &s_t[sOffLagrange];
    PylithScalar nDotLagrange_t = 0.;
    for (PylithInt i = 0; i < spaceDim; ++i) {
        nDotLagrange_t -= n[i] * faultLagrange_t[i];
    }

    f0[fOffp_fault] += porosity * (betaP * (pressureN_t + 2. * pressureFault_t + pressureP_t) / 4. +
                                   betaSigma * nDotLagrange_t) -
                       source;
    if (f0[fOffp_fault] != f0[fOffp_fault]) {
        PetscPrintf(PETSC_COMM_WORLD, "Error in f0p_fault_body_source \n");
    }

} // f0p_fault_body_source


// ----------------------------------------------------------------------
// f1 function for p_fault constraint equation:
// f1p_fault = \kappa_{fx} / (4\mu) \vnabla (p^+ + 2 p^f + p^-)
void
pylith::fekernels::FaultCohesiveKinPoro::f1p_fault(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar f1[]) {
    assert(aOff);
    assert(sOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(f1);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;

    const PylithScalar *vectorPermeability = &a[aOff[i_fault_permeability]];
    // const PylithScalar vectorPermeability[4] = {1.0, 1.0, 0.0, 0.0};
    // const PylithScalar fluidViscosity = 1.0;

    PylithScalar tensorPermeability[spaceDim * spaceDim];
    switch (spaceDim) {
    case 1:
        tensorPermeability[0] = vectorPermeability[0];
        break;
    case 2:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[3];
        tensorPermeability[3] = vectorPermeability[1];
        break;
    case 3:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];
        break;
    default:
        assert(0);
    } // switch

    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    // Pressure_x
    const PylithInt sOffPressureN = sOff[i_pressure];
    const PylithInt sOffPressureP = sOffPressureN + 1;
    const PylithInt sOffPressureFault = sOff[i_fault_pressure];

    const PylithScalar *pressureN_x = &s_x[sOff_x[i_pressure]];
    const PylithScalar *pressureP_x = &s_x[sOff_x[i_pressure] + spaceDim];
    const PylithScalar *pressureFault_x = &s_x[sOff_x[i_fault_pressure]];
    const PylithInt fOffp_fault = 0;

    // Do transformation for gradient
    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {n[1], -n[0]};
#if 1
        const PylithScalar dpdx = pressureFault_x[0]*tanDir[0];
        const PylithScalar dpdy = pressureFault_x[0]*tanDir[1];
        const PylithScalar bodyForceX = 0.0; // porosity*0.0;
        const PylithScalar bodyForceY = 0.0; // porosity*0.0;

        f1[fOffp_fault+0] += (tensorPermeability[0]/fluidViscosity) * (dpdx - bodyForceX) + (tensorPermeability[1]/fluidViscosity) * (dpdy - bodyForceY);
        f1[fOffp_fault+1] += (tensorPermeability[2]/fluidViscosity) * (dpdx - bodyForceX) + (tensorPermeability[3]/fluidViscosity) * (dpdy - bodyForceY);
#else
        f1[fOffp_fault] += tensorPermeability[0] / (4.0 * fluidViscosity) * tanDir[0] * tanDir[0] * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[1] / (4.0 * fluidViscosity) * tanDir[0] * tanDir[1] * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]);
        f1[fOffp_fault + 1] += tensorPermeability[2] / (4.0 * fluidViscosity) * tanDir[1] * tanDir[0] * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[3] / (4.0 * fluidViscosity) * tanDir[1] * tanDir[1] * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]);
#endif
        break;
    } // case 2

    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);
        f1[fOffp_fault] += tensorPermeability[0] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[0] + tanDir2[0] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[1] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[1] + tanDir2[0] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) + tensorPermeability[2] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[2] + tanDir2[0] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]);
        f1[fOffp_fault + 1] += tensorPermeability[3] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[0] + tanDir2[1] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[4] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[1] + tanDir2[1] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) + tensorPermeability[5] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[2] + tanDir2[1] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]);
        f1[fOffp_fault + 2] += tensorPermeability[6] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[0] + tanDir2[2] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[7] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[1] + tanDir2[2] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) + tensorPermeability[8] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[2] + tanDir2[2] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]);
        break;
    } // case 3
    default:
        assert(0);
    }
    for (PylithInt i = 0; i < dim; ++i) {
        if (f1[fOffp_fault + i] != f1[fOffp_fault + i]) {
            PetscPrintf(PETSC_COMM_WORLD, "Error in f1p_fault \n");
        }
    }
} // f1p_fault


// ----------------------------------------------------------------------
// f1 function for p_fault constraint equation:
/* f1p_fault = \tensor{\kappa_{f}} / (4\mu) \tensor{I} \tensor{I'} \cdot \left( \vnabla (p^+ + 2 p^f + p^-)
 * - (f^+ + 2 f^f + f^-) \right)
 *
 * FAULT COHESIVE Face
 */
void
pylith::fekernels::FaultCohesiveKinPoro::f1p_fault_body(const PylithInt dim,
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
                                                        const PylithReal n[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar f1[]) {
    assert(aOff);
    assert(sOff_x);
    assert(s);
    assert(a);
    assert(s_x);
    assert(f1);

    assert(numS >= 5);
    // assert(numA >= 5);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    // Index for solution fields
    const PylithInt i_pressure = 1;
    const PylithInt i_fault_pressure = 4;

    // Index for auxiliary fields
    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;
    const PylithInt i_body_force = numA - 3;

    const PylithScalar *vectorPermeability = &a[aOff[i_fault_permeability]];
    PylithScalar tensorPermeability[spaceDim * spaceDim];
    switch (spaceDim) {
    case 1:
        tensorPermeability[0] = vectorPermeability[0];
        break;
    case 2:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];
        break;
    case 3:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];
        break;
    default:
        assert(0);
    } // switch

    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];
    const PylithScalar *body_force_N = &a[aOff[i_body_force]];
    const PylithScalar *body_force_Fault = &a[aOff[i_body_force]];
    const PylithScalar *body_force_P = &a[aOff[i_body_force]];

    // Pressure_x
    const PylithInt sOffPressureN = sOff[i_pressure];
    const PylithInt sOffPressureP = sOffPressureN + 1;
    const PylithInt sOffPressureFault = sOff[i_fault_pressure];

    const PylithScalar *pressureN_x = &s_x[sOff_x[sOffPressureN]];
    const PylithScalar *pressureP_x = &s_x[sOff_x[sOffPressureP]];
    const PylithScalar *pressureFault_x = &s_x[sOff_x[sOffPressureFault]];
    const PylithInt fOffp_fault = 0;

    // Do transformation for gradient
    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        f1[fOffp_fault] += tensorPermeability[0] / (4.0 * fluidViscosity) * tanDir[0] * tanDir[0] * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[1] / (4.0 * fluidViscosity) * tanDir[0] * tanDir[1] * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]);
        f1[fOffp_fault + 1] += tensorPermeability[2] / (4.0 * fluidViscosity) * tanDir[1] * tanDir[0] * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) + tensorPermeability[3] / (4.0 * fluidViscosity) * tanDir[1] * tanDir[1] * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]);
        break;
    } // case 2

    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);
        f1[fOffp_fault] += tensorPermeability[0] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[0] + tanDir2[0] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) * (body_force_N[0] + 2. * body_force_Fault[0] + body_force_P[0]) + tensorPermeability[1] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[1] + tanDir2[0] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) * (body_force_N[1] + 2. * body_force_Fault[1] + body_force_P[1]) + tensorPermeability[2] / (4.0 * fluidViscosity) * (tanDir1[0] * tanDir1[2] + tanDir2[0] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]) * (body_force_N[2] + 2. * body_force_Fault[2] + body_force_P[2]);
        f1[fOffp_fault + 1] += tensorPermeability[3] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[0] + tanDir2[1] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) * (body_force_N[0] + 2. * body_force_Fault[0] + body_force_P[0]) + tensorPermeability[4] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[1] + tanDir2[1] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) * (body_force_N[1] + 2. * body_force_Fault[1] + body_force_P[1]) + tensorPermeability[5] / (4.0 * fluidViscosity) * (tanDir1[1] * tanDir1[2] + tanDir2[1] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]) * (body_force_N[2] + 2. * body_force_Fault[2] + body_force_P[2]);
        f1[fOffp_fault + 2] += tensorPermeability[6] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[0] + tanDir2[2] * tanDir2[0]) * (pressureN_x[0] + 2. * pressureFault_x[0] + pressureP_x[0]) * (body_force_N[0] + 2. * body_force_Fault[0] + body_force_P[0]) + tensorPermeability[7] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[1] + tanDir2[2] * tanDir2[1]) * (pressureN_x[1] + 2. * pressureFault_x[1] + pressureP_x[1]) * (body_force_N[1] + 2. * body_force_Fault[1] + body_force_P[1]) + tensorPermeability[8] / (4.0 * fluidViscosity) * (tanDir1[2] * tanDir1[2] + tanDir2[2] * tanDir2[2]) * (pressureN_x[2] + 2. * pressureFault_x[2] + pressureP_x[2]) * (body_force_N[2] + 2. * body_force_Fault[2] + body_force_P[2]);
        break;
    } // case 3
    default:
        assert(0);
    }
    for (PylithInt i = 0; i < spaceDim; ++i) {
        if (f1[fOffp_fault + i] != f1[fOffp_fault + i]) {
            PetscPrintf(PETSC_COMM_WORLD, "Error in f1p_fault \n");
        }
    }
} // f1p_fault_body


// ----------------------------------------------------------------------
// Jacobian Functions - JFU
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/* Jf0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0ul_neg(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf0[]) {
    assert(numS >= 5);
    // assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    // const PylithInt gOffP = gOffN + spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[(gOffN + i) * ncols + i] += 1.0;
        // Jf0[(gOffP + i) * ncols + i] += +1.0;
    } // for
} // Jg0ul_neg


// ----------------------------------------------------------------------
/* Jf0 function for integration of the displacement equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0ul_pos(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf0[]) {
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        // Jf0[(gOffN + i) * ncols + i] += -1.0;
        Jf0[(gOffP + i) * ncols + i] += -1.0;
    } // for
} // Jg0ul_pos


// ----------------------------------------------------------------------
// Jacobian Functions - JFP
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Jacobian Functions - JFE
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Jacobian Functions - JFU
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Jacobian Functions - JF\lambda
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/* Jg0 function for integration of the slip constraint equation.
 *
 * Solution fields = [disp(dim), ..., lagrange(dim)]
 * Auxiliary fields = None
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0lu(const PylithInt dim,
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
                                               const PylithReal n[],
                                               const PylithInt numConstants,
                                               const PylithScalar constants[],
                                               PylithScalar Jf0[]) {
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN+spaceDim*spaceDim;
    const PylithInt ncols = spaceDim;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[gOffN+i*ncols+i] += +1.0;
        Jf0[gOffP+i*ncols+i] += -1.0;
    } // for
} // Jg0lu


// ----------------------------------------------------------------------
// Jacobian Functions - JFP_fault
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/* Jf0p_fp function for integration of the displacement equation.
 * [\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2,\phi_f beta^p t_shift /4 - \kappa_{fz}/\mu h^2]
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0p_fp(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar Jf0[]) {
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];

    // Hardcoded test values
    // const PylithScalar porosity = 0.1;
    // const PylithScalar betaP = 1.0;

    const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + 1;

    Jf0[gOffN] += (porosity * betaP * s_tshift) / 4.0;
    Jf0[gOffP] += (porosity * betaP * s_tshift) / 4.0;

} // Jf0p_fp


// ----------------------------------------------------------------------
/* Jf3p_fp function for integration of the displacement equation.
 * [\kappa_{fx} / 4 / \mu \te{I}, \kappa_{fx} / 4 / \mu \te{I}]
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf3p_fp(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar Jf3[]) {
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf3);
    // Index for auxiliary fields
    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;

    const PylithInt spaceDim = dim + 1;

    const PylithScalar *vectorPermeability = &a[aOff[i_fault_permeability]];
    // const PylithScalar vectorPermeability[4] = {1.0, 1.0, 0.0, 0.0};
    // const PylithScalar fluidViscosity = 1.0;

    PylithScalar tensorPermeability[spaceDim * spaceDim];
    switch (spaceDim) {
    case 1:
        tensorPermeability[0] = vectorPermeability[0];
        break;
    case 2:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];
        break;
    case 3:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];
        break;
    default:
        assert(0);
    } // switch

    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt gOffN = 0;
    const PylithInt gOffP = gOffN + spaceDim;

    const PylithInt ncols = 2 * spaceDim;
    const PylithScalar tempConst = fluidViscosity / 4.;
    // Do transformation for gradient
    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        Jf3[0 * ncols + 0] += (tensorPermeability[0] / tempConst) * tanDir[0] * tanDir[0];
        Jf3[0 * ncols + 1] += (tensorPermeability[1] / tempConst) * tanDir[0] * tanDir[1];
        Jf3[1 * ncols + 0] += (tensorPermeability[2] / tempConst) * tanDir[1] * tanDir[0];
        Jf3[1 * ncols + 0] += (tensorPermeability[3] / tempConst) * tanDir[1] * tanDir[1];

        Jf3[0 * ncols + 0] += (tensorPermeability[0] / tempConst) * tanDir[0] * tanDir[0];
        Jf3[0 * ncols + 1] += (tensorPermeability[1] / tempConst) * tanDir[0] * tanDir[1];
        Jf3[1 * ncols + 0] += (tensorPermeability[2] / tempConst) * tanDir[1] * tanDir[0];
        Jf3[1 * ncols + 0] += (tensorPermeability[3] / tempConst) * tanDir[1] * tanDir[1];

        Jf3[0 * ncols + gOffP + 0] += (tensorPermeability[0] / tempConst) * tanDir[0] * tanDir[0];
        Jf3[0 * ncols + gOffP + 1] += (tensorPermeability[1] / tempConst) * tanDir[0] * tanDir[1];
        Jf3[1 * ncols + gOffP + 0] += (tensorPermeability[2] / tempConst) * tanDir[1] * tanDir[0];
        Jf3[1 * ncols + gOffP + 0] += (tensorPermeability[3] / tempConst) * tanDir[1] * tanDir[1];

        Jf3[0 * ncols + gOffP + 0] += (tensorPermeability[0] / tempConst) * tanDir[0] * tanDir[0];
        Jf3[0 * ncols + gOffP + 1] += (tensorPermeability[1] / tempConst) * tanDir[0] * tanDir[1];
        Jf3[1 * ncols + gOffP + 0] += (tensorPermeability[2] / tempConst) * tanDir[1] * tanDir[0];
        Jf3[1 * ncols + gOffP + 0] += (tensorPermeability[3] / tempConst) * tanDir[1] * tanDir[1];
        break;
    } // case 2

    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);
        Jf3[0 * ncols + 0] += (tensorPermeability[0] / tempConst) * (tanDir1[0] * tanDir1[0] + tanDir2[0] * tanDir2[0]);
        Jf3[0 * ncols + 1] += (tensorPermeability[1] / tempConst) * (tanDir1[0] * tanDir1[1] + tanDir2[0] * tanDir2[1]);
        Jf3[0 * ncols + 2] += (tensorPermeability[2] / tempConst) * (tanDir1[0] * tanDir1[2] + tanDir2[0] * tanDir2[2]);
        Jf3[1 * ncols + 0] += (tensorPermeability[3] / tempConst) * (tanDir1[1] * tanDir1[0] + tanDir2[1] * tanDir2[0]);
        Jf3[1 * ncols + 1] += (tensorPermeability[4] / tempConst) * (tanDir1[1] * tanDir1[1] + tanDir2[1] * tanDir2[1]);
        Jf3[1 * ncols + 2] += (tensorPermeability[5] / tempConst) * (tanDir1[1] * tanDir1[2] + tanDir2[1] * tanDir2[2]);
        Jf3[2 * ncols + 0] += (tensorPermeability[6] / tempConst) * (tanDir1[2] * tanDir1[0] + tanDir2[2] * tanDir2[0]);
        Jf3[2 * ncols + 1] += (tensorPermeability[7] / tempConst) * (tanDir1[2] * tanDir1[1] + tanDir2[2] * tanDir2[1]);
        Jf3[2 * ncols + 2] += (tensorPermeability[8] / tempConst) * (tanDir1[2] * tanDir1[2] + tanDir2[2] * tanDir2[2]);

        Jf3[0 * ncols + gOffP + 0] += (tensorPermeability[0] / tempConst) * (tanDir1[0] * tanDir1[0] + tanDir2[0] * tanDir2[0]);
        Jf3[0 * ncols + gOffP + 1] += (tensorPermeability[1] / tempConst) * (tanDir1[0] * tanDir1[1] + tanDir2[0] * tanDir2[1]);
        Jf3[0 * ncols + gOffP + 2] += (tensorPermeability[2] / tempConst) * (tanDir1[0] * tanDir1[2] + tanDir2[0] * tanDir2[2]);
        Jf3[1 * ncols + gOffP + 0] += (tensorPermeability[3] / tempConst) * (tanDir1[1] * tanDir1[0] + tanDir2[1] * tanDir2[0]);
        Jf3[1 * ncols + gOffP + 1] += (tensorPermeability[4] / tempConst) * (tanDir1[1] * tanDir1[1] + tanDir2[1] * tanDir2[1]);
        Jf3[1 * ncols + gOffP + 2] += (tensorPermeability[5] / tempConst) * (tanDir1[1] * tanDir1[2] + tanDir2[1] * tanDir2[2]);
        Jf3[2 * ncols + gOffP + 0] += (tensorPermeability[6] / tempConst) * (tanDir1[2] * tanDir1[0] + tanDir2[2] * tanDir2[0]);
        Jf3[2 * ncols + gOffP + 1] += (tensorPermeability[7] / tempConst) * (tanDir1[2] * tanDir1[1] + tanDir2[2] * tanDir2[1]);
        Jf3[2 * ncols + gOffP + 2] += (tensorPermeability[8] / tempConst) * (tanDir1[2] * tanDir1[2] + tanDir2[2] * tanDir2[2]);
        break;
    } // case 3
    default:
        assert(0);
    }
    /**
     * for (PylithInt i = 0; i < spaceDim; ++i) {
     *  Jf3[i * ncols + gOffN + i] += permeabilityTangential / (4. * fluidViscosity);
     *  Jf3[i * ncols + gOffP + i] += permeabilityTangential / (4. * fluidViscosity);
     * }
     */
} // Jf3p_fp


// ----------------------------------------------------------------------
/* Jf0p_fl function for integration of the displacement equation.
 * s_tshift \phi_f \beta^\sigma \ve{n}
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0p_fl(const PylithInt dim,
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
                                                 const PylithReal n[],
                                                 const PylithInt numConstants,
                                                 const PylithScalar constants[],
                                                 PylithScalar Jf0[]) {
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);

    const PylithInt spaceDim = dim + 1;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_sigma = 3;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaSigma = a[aOff[i_beta_sigma]];

    // Hardcoded test values
    // const PylithScalar porosity = 0.1;
    // const PylithScalar betaSigma = 1.0;

    for (PylithInt i = 0; i < spaceDim; ++i) {
        Jf0[i] += -betaSigma * porosity * s_tshift * n[i];
    }
} // Jf0p_fl


// ----------------------------------------------------------------------
/* Jf0p_fp_f function for integration of the displacement equation:
 * 2 \kappa_{fz} / (\mu h^2) + 2 \phi_f \beta^p t_shift
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf0p_fp_f(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf0[]) {
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf0);
    assert(sOff);
    assert(aOff);
    assert(n);
    const PylithInt gOff = 0;

    // Index for auxiliary fields
    const PylithInt i_porosity = 1;
    const PylithInt i_beta_p = 2;

    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar betaP = a[aOff[i_beta_p]];

    // Hardcoded test values
    // const PylithScalar porosity = 0.1;
    // const PylithScalar betaP = 1.0;

    Jf0[gOff] += (2.0 * porosity * betaP * s_tshift) / 4.0;

} // Jf0p_fp_f


// ----------------------------------------------------------------------
/* Jf3p_fp_f function for integration of the displacement equation:
 * \kappa_{fx} / (2 \mu) \te{I}
 * Solution fields = [disp(dim), ..., lagrange(dim), fault_pressure(1)]
 * Auxiliary fields
 */
void
pylith::fekernels::FaultCohesiveKinPoro::Jf3p_fp_f(const PylithInt dim,
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
                                                   const PylithReal n[],
                                                   const PylithInt numConstants,
                                                   const PylithScalar constants[],
                                                   PylithScalar Jf3[]) {
    // Check data fields
    assert(numS >= 5);
    assert(numA >= 5);
    assert(Jf3);
    assert(sOff);
    assert(aOff);
    assert(n);
    const PylithInt gOff = 0;

    const PylithInt spaceDim = dim + 1;

    // Index for auxiliary fields
    const PylithInt i_fault_permeability = 4;
    const PylithInt i_fluid_viscosity = 5;

    const PylithScalar *vectorPermeability = &a[aOff[i_fault_permeability]];

    // Hardcoded test values
    // const PylithScalar vectorPermeability[4] = {1.0, 1.0, 0.0, 0.0};
    // const PylithScalar fluidViscosity = 1.0;

    PylithScalar tensorPermeability[spaceDim * spaceDim];
    switch (spaceDim) {
    case 1:
        tensorPermeability[0] = vectorPermeability[0];
        break;
    case 2:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];
        break;
    case 3:
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];
        break;
    default:
        assert(0);
    } // switch

    const PylithScalar fluidViscosity = a[aOff[i_fluid_viscosity]];

    const PylithInt ncols = spaceDim;
    const PylithScalar tempConst = 2.0 / (4.0 * fluidViscosity);
    // Do transformation for gradient
    switch (spaceDim) {
    case 2:
    {
        const PylithInt _spaceDim = 2;
        const PylithScalar tanDir[2] = {-n[1], n[0]};
        Jf3[0 * ncols + 0] += tensorPermeability[0] * tempConst * tanDir[0] * tanDir[0];
        Jf3[0 * ncols + 1] += tensorPermeability[1] * tempConst * tanDir[0] * tanDir[1];
        Jf3[1 * ncols + 0] += tensorPermeability[2] * tempConst * tanDir[1] * tanDir[0];
        Jf3[1 * ncols + 0] += tensorPermeability[3] * tempConst * tanDir[1] * tanDir[1];
        break;
    } // case 2

    case 3:
    {
        const PylithInt _spaceDim = 3;
        const PylithScalar *refDir1 = &constants[0];
        const PylithScalar *refDir2 = &constants[3];
        PylithScalar tanDir1[3], tanDir2[3];
        pylith::fekernels::_FaultCohesiveKinPoro::tangential_directions(_spaceDim, refDir1, refDir2, n, tanDir1, tanDir2);
        Jf3[0 * ncols + 0] += tensorPermeability[0] * tempConst * (tanDir1[0] * tanDir1[0] + tanDir2[0] * tanDir2[0]);
        Jf3[0 * ncols + 1] += tensorPermeability[1] * tempConst * (tanDir1[0] * tanDir1[1] + tanDir2[0] * tanDir2[1]);
        Jf3[0 * ncols + 2] += tensorPermeability[2] * tempConst * (tanDir1[0] * tanDir1[2] + tanDir2[0] * tanDir2[2]);
        Jf3[1 * ncols + 0] += tensorPermeability[3] * tempConst * (tanDir1[1] * tanDir1[0] + tanDir2[1] * tanDir2[0]);
        Jf3[1 * ncols + 1] += tensorPermeability[4] * tempConst * (tanDir1[1] * tanDir1[1] + tanDir2[1] * tanDir2[1]);
        Jf3[1 * ncols + 2] += tensorPermeability[5] * tempConst * (tanDir1[1] * tanDir1[2] + tanDir2[1] * tanDir2[2]);
        Jf3[2 * ncols + 0] += tensorPermeability[6] * tempConst * (tanDir1[2] * tanDir1[0] + tanDir2[2] * tanDir2[0]);
        Jf3[2 * ncols + 1] += tensorPermeability[7] * tempConst * (tanDir1[2] * tanDir1[1] + tanDir2[2] * tanDir2[1]);
        Jf3[2 * ncols + 2] += tensorPermeability[8] * tempConst * (tanDir1[2] * tanDir1[2] + tanDir2[2] * tanDir2[2]);
        break;
    } // case 3

    default:
        assert(0);
    }

} // Jf3p_fp_f


// End of file
