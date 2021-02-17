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

#include "pylith/fekernels/PointSource.hh"

#include <cmath> // USES exp()
#include <cassert> // USES assert()
#include <iostream> // debugging.


// ----------------------------------------------------------------------
// g0 function for displacement point source, ricker source time function
void
pylith::fekernels::PointSource::g0u_ricker(const PylithInt dim,
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
    assert(sOff);
    assert(s_x);
    assert(g0);

    const PylithInt i_disp = 1;
    const PylithScalar* disp_x = &s_x[sOff[i_vel]];
    const PylithScalar time_shift = constants[9];
    const PylithScalar f_0 = constants[10];
    const PylithScalar t_adj = t-time_shift;

    PylithScalar momentTensor[9];

    momentTensor[0] = constants[0]; // Mrr / Mxx
    momentTensor[1] = constants[3]; // Mrt / Mxy
    momentTensor[2] = constants[4]; // Mrp / Mxz
    momentTensor[3] = constants[3]; // Mrt / Mxy
    momentTensor[4] = constants[1]; // Mtt / Myy
    momentTensor[5] = constants[5]; // Mtp / Myz
    momentTensor[6] = constants[4]; // Mrp / Mxz
    momentTensor[7] = constants[5]; // Mtp / Myz
    momentTensor[8] = constants[2]; // Mpp / Mzz

    const PylithScalar rickerAmplitude = (1.0 - 2.0*M_PI*M_PI*f_0*f_0*t_adj*t_adj)*exp(-1.0*M_PI*M_PI*f_0*f_0*t_adj*t_adj);

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; ++j) {
            g0[i] += -1.0*rickerAmplitude*momentTensor[i*dim+j]*disp_x[j];
        } // for
    } // for
} // g0u_ricker

// End of file
