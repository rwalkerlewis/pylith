/*
 * =================================================================================================
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
 * Point force kernels for elastodynamics.
 *
 * The point force is implemented using a regularized delta function (Gaussian).
 * The force contribution at point x from a moment tensor source at xs is:
 *
 * f_i(x, t) = M_ij * dphi/dx_j * S(t) * magnitude
 *
 * where:
 *   M_ij is the moment tensor
 *   phi is the regularized delta function: phi(r) = exp(-r^2 / (2*sigma^2)) / (2*pi*sigma^2)^(dim/2)
 *   dphi/dx_j is the spatial derivative
 *   S(t) is the source time function (Ricker wavelet)
 *   magnitude is the seismic moment
 *
 * Kernel interface:
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
 * @param[out] f0 [dim].
 */

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cmath> // USES exp, sin, cos, sqrt

// ------------------------------------------------------------------------------------------------
/// Kernels for point force in 2D.
class pylith::fekernels::PointForce2D {
public:

    // --------------------------------------------------------------------------------------------
    /** Ricker wavelet source time function.
     *
     * S(t) = (1 - 2 * pi^2 * f^2 * tau^2) * exp(-pi^2 * f^2 * tau^2)
     *
     * where tau = t - t0 - delay, f is dominant frequency, t0 is origin time.
     */
    static inline
    PylithReal rickerWavelet(const PylithReal t,
                             const PylithReal originTime,
                             const PylithReal frequency,
                             const PylithReal delay) {
        const PylithReal pi = M_PI;
        const PylithReal tau = t - originTime - delay;
        const PylithReal pi2_f2_tau2 = (pi * frequency * tau) * (pi * frequency * tau);
        return (1.0 - 2.0 * pi2_f2_tau2) * exp(-pi2_f2_tau2);
    } // rickerWavelet

    // --------------------------------------------------------------------------------------------
    /** Regularized delta function (2D Gaussian).
     *
     * phi(r) = exp(-r^2 / (2*sigma^2)) / (2*pi*sigma^2)
     */
    static inline
    PylithReal deltaFunction2D(const PylithReal dist,
                               const PylithReal sigma) {
        const PylithReal pi = M_PI;
        const PylithReal norm = 1.0 / (2.0 * pi * sigma * sigma);
        return norm * exp(-dist * dist / (2.0 * sigma * sigma));
    } // deltaFunction2D

    // --------------------------------------------------------------------------------------------
    /** Gradient of regularized delta function (2D).
     *
     * dphi/dr = -r / sigma^2 * phi(r)
     */
    static inline
    void deltaGradient2D(const PylithReal* r,
                         const PylithReal dist,
                         const PylithReal sigma,
                         PylithReal* grad) {
        const PylithReal phi = deltaFunction2D(dist, sigma);
        const PylithReal factor = -1.0 / (sigma * sigma) * phi;
        grad[0] = r[0] * factor;
        grad[1] = r[1] * factor;
    } // deltaGradient2D

    // --------------------------------------------------------------------------------------------
    /** g0 function for point force contribution to velocity equation.
     *
     * ISA PetscPointFn*
     *
     * Auxiliary fields: [source_location(2), moment_tensor(3), magnitude, origin_time, frequency, delay]
     * Constants: [origin_time, frequency, delay, magnitude]
     */
    static inline
    void g0v_pointforce(const PylithInt dim,
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
        const PylithInt _dim = 2;
        const PylithInt _numA = 6;
        const PylithInt _numConstants = 4;

        assert(_dim == dim);
        assert(_numA <= numA);
        assert(_numConstants <= numConstants);
        assert(aOff);
        assert(a);
        assert(constants);
        assert(g0);

        // Auxiliary field indices
        const PylithInt i_sourceLocation = 0;
        const PylithInt i_momentTensor = 1;
        const PylithInt i_magnitude = 2;
        const PylithInt i_originTime = 3;
        const PylithInt i_frequency = 4;
        const PylithInt i_delay = 5;

        // Get auxiliary field values
        const PylithScalar* sourceLocation = &a[aOff[i_sourceLocation]];
        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]]; // [Mxx, Myy, Mxy]
        const PylithScalar magnitude = a[aOff[i_magnitude]];
        const PylithScalar originTime = a[aOff[i_originTime]];
        const PylithScalar frequency = a[aOff[i_frequency]];
        const PylithScalar delay = a[aOff[i_delay]];

        // Compute distance from source
        PylithReal r[2];
        r[0] = x[0] - sourceLocation[0];
        r[1] = x[1] - sourceLocation[1];
        const PylithReal dist = sqrt(r[0]*r[0] + r[1]*r[1]);

        // Characteristic length (regularization parameter)
        // Use a fraction of the domain size or element size
        const PylithReal sigma = 0.01; // Should be set based on mesh resolution

        // Source time function (Ricker wavelet)
        const PylithReal stf = rickerWavelet(t, originTime, frequency, delay);

        // Gradient of delta function
        PylithReal deltaGrad[2] = {0.0, 0.0};
        if (dist > 1.0e-12 * sigma) {
            deltaGradient2D(r, dist, sigma, deltaGrad);
        }

        // Moment tensor in matrix form:
        // M = [Mxx  Mxy]
        //     [Mxy  Myy]
        const PylithReal Mxx = momentTensor[0];
        const PylithReal Myy = momentTensor[1];
        const PylithReal Mxy = momentTensor[2];

        // Compute equivalent force: f_i = -M_ij * dphi/dx_j * S(t) * magnitude
        // f_x = -(Mxx * dphi/dx + Mxy * dphi/dy) * S(t) * magnitude
        // f_y = -(Mxy * dphi/dx + Myy * dphi/dy) * S(t) * magnitude
        g0[0] += -(Mxx * deltaGrad[0] + Mxy * deltaGrad[1]) * stf * magnitude;
        g0[1] += -(Mxy * deltaGrad[0] + Myy * deltaGrad[1]) * stf * magnitude;
    } // g0v_pointforce

}; // PointForce2D


// ------------------------------------------------------------------------------------------------
/// Kernels for point force in 3D.
class pylith::fekernels::PointForce3D {
public:

    // --------------------------------------------------------------------------------------------
    /** Ricker wavelet source time function. */
    static inline
    PylithReal rickerWavelet(const PylithReal t,
                             const PylithReal originTime,
                             const PylithReal frequency,
                             const PylithReal delay) {
        return PointForce2D::rickerWavelet(t, originTime, frequency, delay);
    } // rickerWavelet

    // --------------------------------------------------------------------------------------------
    /** Regularized delta function (3D Gaussian).
     *
     * phi(r) = exp(-r^2 / (2*sigma^2)) / (2*pi*sigma^2)^(3/2)
     */
    static inline
    PylithReal deltaFunction3D(const PylithReal dist,
                               const PylithReal sigma) {
        const PylithReal pi = M_PI;
        const PylithReal norm = 1.0 / pow(2.0 * pi * sigma * sigma, 1.5);
        return norm * exp(-dist * dist / (2.0 * sigma * sigma));
    } // deltaFunction3D

    // --------------------------------------------------------------------------------------------
    /** Gradient of regularized delta function (3D). */
    static inline
    void deltaGradient3D(const PylithReal* r,
                         const PylithReal dist,
                         const PylithReal sigma,
                         PylithReal* grad) {
        const PylithReal phi = deltaFunction3D(dist, sigma);
        const PylithReal factor = -1.0 / (sigma * sigma) * phi;
        grad[0] = r[0] * factor;
        grad[1] = r[1] * factor;
        grad[2] = r[2] * factor;
    } // deltaGradient3D

    // --------------------------------------------------------------------------------------------
    /** g0 function for point force contribution to velocity equation.
     *
     * ISA PetscPointFn*
     *
     * Auxiliary fields: [source_location(3), moment_tensor(6), magnitude, origin_time, frequency, delay]
     */
    static inline
    void g0v_pointforce(const PylithInt dim,
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
        const PylithInt _dim = 3;
        const PylithInt _numA = 6;

        assert(_dim == dim);
        assert(_numA <= numA);
        assert(aOff);
        assert(a);
        assert(g0);

        // Auxiliary field indices
        const PylithInt i_sourceLocation = 0;
        const PylithInt i_momentTensor = 1;
        const PylithInt i_magnitude = 2;
        const PylithInt i_originTime = 3;
        const PylithInt i_frequency = 4;
        const PylithInt i_delay = 5;

        // Get auxiliary field values
        const PylithScalar* sourceLocation = &a[aOff[i_sourceLocation]];
        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]]; // [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
        const PylithScalar magnitude = a[aOff[i_magnitude]];
        const PylithScalar originTime = a[aOff[i_originTime]];
        const PylithScalar frequency = a[aOff[i_frequency]];
        const PylithScalar delay = a[aOff[i_delay]];

        // Compute distance from source
        PylithReal r[3];
        r[0] = x[0] - sourceLocation[0];
        r[1] = x[1] - sourceLocation[1];
        r[2] = x[2] - sourceLocation[2];
        const PylithReal dist = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

        // Characteristic length (regularization parameter)
        const PylithReal sigma = 0.01;

        // Source time function (Ricker wavelet)
        const PylithReal stf = rickerWavelet(t, originTime, frequency, delay);

        // Gradient of delta function
        PylithReal deltaGrad[3] = {0.0, 0.0, 0.0};
        if (dist > 1.0e-12 * sigma) {
            deltaGradient3D(r, dist, sigma, deltaGrad);
        }

        // Moment tensor components: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
        const PylithReal Mxx = momentTensor[0];
        const PylithReal Myy = momentTensor[1];
        const PylithReal Mzz = momentTensor[2];
        const PylithReal Mxy = momentTensor[3];
        const PylithReal Mxz = momentTensor[4];
        const PylithReal Myz = momentTensor[5];

        // Compute equivalent force: f_i = -M_ij * dphi/dx_j * S(t) * magnitude
        g0[0] += -(Mxx * deltaGrad[0] + Mxy * deltaGrad[1] + Mxz * deltaGrad[2]) * stf * magnitude;
        g0[1] += -(Mxy * deltaGrad[0] + Myy * deltaGrad[1] + Myz * deltaGrad[2]) * stf * magnitude;
        g0[2] += -(Mxz * deltaGrad[0] + Myz * deltaGrad[1] + Mzz * deltaGrad[2]) * stf * magnitude;
    } // g0v_pointforce

}; // PointForce3D

// End of file
