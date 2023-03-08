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

/** @file libsrc/fekernels/IsotropicLinearPoroelasticity.hh
 *
 * Kernels for linear poroelasticity plane strain.
 *
 * Solution fields: [disp(dim), pressure(1),trace_strain(1) ] (QS)
 * OR
 * Solution fields: [disp(dim), pressure(1),velocity(dim) ] (QS)
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 ***** Required fields(govening equations) + option fields + required fields (rheology)
 * - 0: solid_density(1)
 * - 1: fluid_density(1)
 * - 2: fluid_viscosity(1)
 * - 3: porosity(1)
 *
 ** Optional fields
 * - +1: gravity_field (dim, optional)
 * - +1: body_force(dim,optional)
 * - +1: source_density(1,optional)
 * - +1: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - +1: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 ** Rheological fields
 * - numA - 5: addShearModulus(1)
 * - numA - 4: addDrainedBulkModulus(1)
 * - numA - 3: addBiotCoefficient(1)
 * - numA - 2: addIsotropicPermeability(1)
 * - numA - 1: addFluidBulkModulus(1)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ======================================================================
 *
 * Kernel interface.
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

#if !defined(pylith_fekernels_isotropiclinearporoelasticity_hh)
#define pylith_fekernels_isotropiclinearporoelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity.
class pylith::fekernels::IsotropicLinearPoroelasticity {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal pressure;
        PylithReal solidDensity;
        PylithReal fluidDensity;
        PylithReal fluidViscosity;
        PylithReal porosity;
        PylithReal shearModulus;
        PylithReal drainedBulkModulus;
        PylithReal biotCoefficient;
        PylithReal fluidBulkModulus;
        pylith::fekernels::Tensor refStress;
        pylith::fekernels::Tensor refStrain;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

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

        // Incoming solution fields.
        const PylithInt i_pressure = 1;
        // Incoming auxiliary fields.
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_porosity = 2;
        const PylithInt i_fluidViscosity = 3;

        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_fluidBulkModulus = numA - 1;

        assert(numA >= 3); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(s);
        assert(sOff);
        assert(sOff[i_pressure]);

        context->pressure = s[sOff[i_pressure]];assert(context->pressure > 0.0);

        context->solidDensity = a[aOff[i_solidDensity]];assert(context->solidDensity > 0.0);
        context->fluidDensity = a[aOff[i_fluidDensity]];assert(context->fluidDensity > 0.0);
        context->porosity = a[aOff[i_porosity]];assert(context->porosity >= 0.0);
        context->fluidViscosity = a[aOff[i_fluidViscosity]];assert(context->fluidViscosity = 0.0);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->drainedBulkModulus > 0.0);
        context->biotCoefficient = a[aOff[i_biotCoefficient]];assert(context->biotCoefficient > 0.0);
        context->fluidBulkModulus = a[aOff[i_fluidBulkModulus]];assert(context->fluidBulkModulus > 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext_refState(Context* context,
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
        const PylithInt i_pressure = 1;
        // Incoming auxiliary fields.
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_porosity = 2;
        const PylithInt i_fluidViscosity = 3;

        const PylithInt i_refStress = numA-7;
        const PylithInt i_refStrain = numA-6;
        const PylithInt i_shearModulus = numA-5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_fluidBulkModulus = numA - 1;

        assert(numA >= 5); // also have density
        assert(a);
        assert(aOff);
        assert(aOff[i_refStress] >= 0);
        assert(aOff[i_refStrain] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(s);
        assert(sOff);
        assert(sOff[i_pressure]);

        context->pressure = s[sOff[i_pressure]];assert(context->pressure > 0.0);

        context->solidDensity = a[aOff[i_solidDensity]];assert(context->solidDensity > 0.0);
        context->fluidDensity = a[aOff[i_fluidDensity]];assert(context->fluidDensity > 0.0);
        context->porosity = a[aOff[i_porosity]];assert(context->porosity >= 0.0);
        context->fluidViscosity = a[aOff[i_fluidViscosity]];assert(context->fluidViscosity = 0.0);

        context->shearModulus = a[aOff[i_shearModulus]];assert(context->shearModulus > 0.0);
        context->drainedBulkModulus = a[aOff[i_drainedBulkModulus]];assert(context->drainedBulkModulus > 0.0);
        context->biotCoefficient = a[aOff[i_biotCoefficient]];assert(context->biotCoefficient > 0.0);
        context->fluidBulkModulus = a[aOff[i_fluidBulkModulus]];assert(context->fluidBulkModulus > 0.0);

        tensorOps.fromVector(&a[aOff[i_refStress]], &context->refStress);
        tensorOps.fromVector(&a[aOff[i_refStrain]], &context->refStrain);
    } // createContext

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress for WITHOUT a reference stress and strain.
     *
     * ISA Poroelasticity::stressFn
     *
     * @param[in] rheologyContext IsotropicLinearPoroelasticity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress(void* rheologyContext,
                      const pylith::fekernels::Tensor& strain,
                      const pylith::fekernels::TensorOps& tensorOps,
                      pylith::fekernels::Tensor* stress) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(stress);

        meanStress(context->pressure, context->drainedBulkModulus, context->biotCoefficient, strain, stress);
        deviatoricStress(context->shearModulus, strain, stress);
    } // cauchyStress

    // --------------------------------------------------------------------------------------------
    /** Helper function for calculating Cauchy stress WITH reference stress/strain.
     *
     * @param[in] rheologyContext IsotropicLinearPoroelastcity context.
     * @param[in] strain Strain tensor.
     * @param[in] tensorOps Tensor operations.
     * @param[out] stress Stress tensor.
     */
    static inline
    void cauchyStress_refState(void* rheologyContext,
                               const pylith::fekernels::Tensor& strain,
                               const pylith::fekernels::TensorOps& tensorOps,
                               pylith::fekernels::Tensor* stress) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(stress);

        const pylith::fekernels::Tensor& refStress = context->refStress;
        const pylith::fekernels::Tensor& refStrain = context->refStrain;
        meanStress_refState(context->pressure, context->drainedBulkModulus, context->biotCoefficient, refStress, refStrain, strain, stress);
        deviatoricStress_refState(context->shearModulus, refStress, refStrain, strain, stress);
    } // cauchyStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate mean (volumetric) stress WITHOUT reference stress and reference strain.
     */
    static inline
    void meanStress(const PylithReal pressure,
                    const PylithReal drainedBulkModulus,
                    const PylithReal biotCoefficient,
                    const pylith::fekernels::Tensor& strain,
                    pylith::fekernels::Tensor* stress) {
        assert(drainedBulkModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal meanStress = drainedBulkModulus * strainTrace - biotCoefficient*pressure;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITHOUT reference stress and strain.
     */
    static inline
    void deviatoricStress(const PylithReal shearModulus,
                          const pylith::fekernels::Tensor& strain,
                          pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;

        stress->xx += 2.0*shearModulus*strain.xx + traceTerm;
        stress->yy += 2.0*shearModulus*strain.yy + traceTerm;
        stress->zz += 2.0*shearModulus*strain.zz + traceTerm;
        stress->xy += 2.0*shearModulus*strain.xy;
        stress->yz += 2.0*shearModulus*strain.yz;
        stress->xz += 2.0*shearModulus*strain.xz;
    } // deviatoricStress

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress WITH reference stress and reference strain.
     */
    static inline
    void meanStress_refState(const PylithReal pressure,
                             const PylithReal drainedBulkModulus,
                             const PylithReal biotCoefficient,
                             const pylith::fekernels::Tensor& refStress,
                             const pylith::fekernels::Tensor& refStrain,
                             const pylith::fekernels::Tensor& strain,
                             pylith::fekernels::Tensor* stress) {
        // Incoming auxiliary fields.
        assert(drainedBulkModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;

        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal meanStress = meanRefStress + drainedBulkModulus * (strainTrace - refStrainTrace) - biotCoefficient*pressure;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress_refState

    // --------------------------------------------------------------------------------------------
    /** Calculate deviatoric stress WITH reference stress and strain.
     */
    static inline
    void deviatoricStress_refState(const PylithReal shearModulus,
                                   const pylith::fekernels::Tensor& refStress,
                                   const pylith::fekernels::Tensor& refStrain,
                                   const pylith::fekernels::Tensor& strain,
                                   pylith::fekernels::Tensor* stress) {
        assert(shearModulus > 0.0);
        assert(stress);

        const PylithReal strainTrace = strain.xx + strain.yy + strain.zz;
        const PylithReal refStrainTrace = refStrain.xx + refStrain.yy + refStrain.zz;
        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refStrainTrace);

        stress->xx += refStress.xx - meanRefStress + 2.0*shearModulus*(strain.xx-refStrain.xx) + traceTerm;
        stress->yy += refStress.yy - meanRefStress + 2.0*shearModulus*(strain.yy-refStrain.yy) + traceTerm;
        stress->zz += refStress.zz - meanRefStress + 2.0*shearModulus*(strain.zz-refStrain.zz) + traceTerm;
        stress->xy += refStress.xy + 2.0*shearModulus*(strain.xy - refStrain.xy);
        stress->yz += refStress.yz + 2.0*shearModulus*(strain.yz - refStrain.yz);
        stress->xz += refStress.xz + 2.0*shearModulus*(strain.xz - refStrain.xz);
    } // deviatoricStress_refState

}; // IsotropicLinearPoroelasticity

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity in Plane Strain.
class pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // ================================= MMS =======================================

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and linear time MMS.
    static inline
    void f0_mms_ql_u(const PylithInt dim,
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
        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        f0[0] -= (2.0 * shearModulus - biotCoefficient * t);
        f0[1] -= (2.0 * lambda + 4.0 * shearModulus - biotCoefficient * t);
    } // f0_quadratic_linear_u

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and linear time MMS.
    static inline
    void f0_mms_ql_p(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff_x);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        PylithScalar sum = 0.0;
        sum += x[0];
        sum += x[1];

        f0[0] -= (sum / biotModulus);
    } // f0_quadratic_linear_p

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and trigonometric time MMS.
    static inline
    void f0_mms_qt_u(const PylithInt dim,
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
        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        f0[0] -= (2.0 * shearModulus - biotCoefficient * PetscCosReal(t));
        f0[1] -= (2.0 * lambda + 4.0 * shearModulus - biotCoefficient * PetscCosReal(t));
    } // f0_quadratic_trig_u

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and trigonometric time MMS.
    static inline
    void f0_mms_qt_p(const PylithInt dim,
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
        // Incoming re-packed solution field.

        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        PylithScalar sum = 0.0;

        sum += x[0];
        sum += x[1];

        f0[0] += PetscSinReal(t) * sum / biotModulus;
    } // f0_quadratic_trig_p

    // ----------------------------------------------------------------------
    // f0u function for trigonometric space and linear time MMS.
    static inline
    void f0_mms_tl_u(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        // f0[0] += PetscSqr(2.0 * PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[0]) * (2.0 * shearModulus + lambda) + 2.0
        // *
        //          (shearModulus + lambda) - 2.0 * PETSC_PI * biotCoefficient * PetscSinReal(2.0 * PETSC_PI * x[0]) *
        // t;
        // f0[1] += PetscSqr(2.0 * PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[1]) * (2.0 * shearModulus + lambda) - 2.0
        // *
        //          PETSC_PI * biotCoefficient * PetscSinReal(2.0 * PETSC_PI * x[1]) * t;

        for (PylithInt d = 0; d < _dim - 1; ++d) {
            f0[d] += PetscSqr(2. * PETSC_PI) * PetscSinReal(2. * PETSC_PI * x[d]) * (2. * shearModulus + lambda) + 2.0 * (shearModulus + lambda) - 2. * PETSC_PI * biotCoefficient * PetscSinReal(2. * PETSC_PI * x[d]) * t;
        }
        f0[_dim - 1] += PetscSqr(2. * PETSC_PI) * PetscSinReal(2. * PETSC_PI * x[_dim - 1]) * (2. * shearModulus + lambda) - 2. * PETSC_PI * biotCoefficient * PetscSinReal(2. * PETSC_PI * x[_dim - 1]) * t;
    } // f0_trig_linear_u

    // ----------------------------------------------------------------------
    // f0p function for trigonometric space and linear time MMS.
    static inline
    void f0_mms_tl_p(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;
        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];
        const PylithScalar kappa = a[aOff[i_isotropicPermeability]] / a[aOff[i_fluidViscosity]];

        PylithScalar sum = 0.0;

        // sum += PetscCosReal(2.0 * PETSC_PI * x[0]) + PetscCosReal(2.0 * PETSC_PI * x[1]);
        for (PylithInt d = 0; d < _dim; ++d) {
            sum += PetscCosReal(2. * PETSC_PI * x[d]);
        }

        f0[0] -= sum / biotModulus - 4.0 * PetscSqr(PETSC_PI) * kappa * sum * t;
    } // f0_quadratic_trig_p

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_explicit(const PylithInt dim,
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
        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += pressure_t / biotModulus;
    } // f0p_explicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        // Incoming re-packed auxiliary field.
        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += s_t ? (biotCoefficient * trace_strain_t) : 0.0;
        f0[0] += s_t ? (pressure_t / biotModulus) : 0.0;
    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 4;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static
    void f0p_implicit_source_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 5;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 5;

        // Incoming re-packed auxiliary field.
        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar source = a[aOff[i_source]];

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 6;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav_body

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
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

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(sOff);
        assert(sOff[i_displacement] >= 0);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_displacement] >= 0);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(f1);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                f1[c * _dim + d] -= shearModulus * (displacement_x[c * _dim + d] + displacement_x[d * _dim + c]);
            } // for
            f1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * trace_strain;
            f1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // f1u

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void f1u_refstate(const PylithInt dim,
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

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_rstress = numA - 7;
        const PylithInt i_rstrain = numA - 6;
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar *refStressVector = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
        const PylithScalar *refStrainVector = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

        // const PylithScalar ref_trace_strain = refStrainVector[0] + refStrainVector[1] + refStrainVector[2];
        // const PylithScalar mean_ref_stress = (refStressVector[0] + refStressVector[1] + refStressVector[2]) / 3.0;

        // Convert reference vectors to reference tensors (2D)
        PylithScalar refStressTensor[_dim * _dim];
        PylithScalar refStrainTensor[_dim * _dim];

        refStressTensor[0] = refStressVector[0];
        refStressTensor[1] = refStressVector[3];
        refStressTensor[2] = refStressVector[3];
        refStressTensor[3] = refStressVector[1];

        refStrainTensor[0] = refStrainVector[0];
        refStrainTensor[1] = refStrainVector[3];
        refStrainTensor[2] = refStrainVector[3];
        refStrainTensor[3] = refStrainVector[1];

        // for (PylithInt c = 0; c < _dim; ++c) {
        //     for (PylithInt d = 0; d < _dim; ++d) {
        //         f1[c * _dim + d] -= refStressTensor[c * _dim + d] + 2.0 * shearModulus * ((displacement_x[c * _dim +
        // d] +
        // displacement_x[d * _dim + c]) / 2.0 - refStrainTensor[c * _dim + d]);
        //     } // for
        //     f1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * (trace_strain -
        // ref_trace_strain);
        //     // Biot Effective Stress Pressure Correction
        //     f1[c * _dim + c] += biotCoefficient * pressure;
        // } // for

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                f1[c * _dim + d] -= refStressTensor[c * _dim + d] + 2.0 * shearModulus * ( (displacement_x[c * _dim + d] + displacement_x[d * _dim + c]) / 2.0 - refStrainTensor[c * _dim + d] );
            } // for
            f1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * trace_strain;
            f1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        f1[0] += ((isotropicPermeability / fluidViscosity) * pressure_x[0]);
        f1[1] += ((isotropicPermeability / fluidViscosity) * pressure_x[1]);
    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j]);
            } // for
        } // for
    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        const PylithInt i_body_force = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - body_force[d]);
        } // for

    } // f1p_body

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_body_force = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j]);
            } // for
        } // for

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(aOff[i_gravityField] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
        } // for

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_gravityField] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
            } // for
        } // for

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluid_density = 1;
        const PylithInt i_fluid_viscosity = 2;

        const PylithInt i_body_force = 4;
        const PylithInt i_gravity_field = 5;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropic_permeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluid_density] >= 0);
        assert(aOff[i_fluid_viscosity] >= 0);
        assert(aOff[i_isotropic_permeability] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_gravity_field] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluid_density = a[aOff[i_fluid_density]];
        const PylithScalar fluid_viscosity = a[aOff[i_fluid_viscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

        const PylithScalar isotropic_permeability = a[aOff[i_isotropic_permeability]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropic_permeability / fluid_viscosity) * (pressure_x[d] - body_force[d] - fluid_density * gravity_field[d]);
        } // for

    } // f1p_body_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;

        const PylithInt i_body_force = 4;
        const PylithInt i_gravity_field = 5;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_gravity_field] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j] - fluidDensity * gravity_field[j]);
            } // for
        } // for

    } // f1p_body_gravity_tensor_permeability

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 2;

        // Incoming solution field.

        // Incoming auxiliary fields.

        // Isotropic Linear Poroelasticity
        const PylithInt i_shearModulus = numA - 5;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            }
        }

    } // Jf3uu

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf2);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] += biotCoefficient;
        } // for
    } // Jf2up

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(Jf2);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] -= lambda;
        } // for
    } // Jf2ue

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 2;

        // index of Incoming auxiliary fields.
        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(Jf3);

        const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf3[d * _dim + d] += isotropicPermeablity / fluidViscosity;
        }
    } // Jf3pp

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 2;

        // index of Incoming auxiliary fields.
        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(Jf3);

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for
    } // Jf3pp_tensorPermeability

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(Jf0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        Jf0[0] += utshift / biotModulus;
    } // Jf0pp

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 2;

        // Incoming re-packed auxiliary field.
        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf0);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        Jf0[0] += utshift * biotCoefficient;
    } // Jf0pe

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0ppdot(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(Jf0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        Jf0[0] += 1.0 / biotModulus;
    } // Jf0ppdot

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pedot(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming re-packed auxiliary field.
        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf0);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        Jf0[0] += biotCoefficient;
    } // Jf0pedot

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_implicit

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
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
        const PylithInt _dim = 2;

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 4;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 5;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_body

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 5;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 6;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav_body

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        for (PylithInt d = 0; d < _dim; ++d) {
            g1[d] -= (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
        } // for

    } // g1p_gravity

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                g1[i] -= (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
            } // for
        } // for

    } // g1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        for (PylithInt d = 0; d < _dim; ++d) {
            g1[d] -= (isotropicPermeability / fluidViscosity) * pressure_x[d];
        } // for
    } // g1p

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        PylithScalar tensorPermeability[4];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[1];
        tensorPermeability[3] = vectorPermeability[3];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                g1[i] -= (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j]);
            } // for
        } // for
    } // g1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void g1v(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain += displacement_x[d * _dim + d];
        }

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                g1[c * dim + d] -= shearModulus * (displacement_x[c * _dim + d] + displacement_x[d * _dim + c]);
            } // for
            g1[c * dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * trace_strain;
            g1[c * dim + c] += biotCoefficient * pressure;
        } // for
    } // g1v

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void g1v_refstate(const PylithInt dim,
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
        const PylithInt _dim = 2;

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain += displacement_x[d * _dim + d];
        }

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_rstress = numA - 7;
        const PylithInt i_rstrain = numA - 6;
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar *refStressVector = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
        const PylithScalar *refStrainVector = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

        const PylithScalar ref_trace_strain = refStrainVector[0] + refStrainVector[1] + refStrainVector[2];
        // const PylithScalar mean_ref_stress = (refStressVector[0] + refStressVector[1] + refStressVector[2]) / 3.0;

        // Convert reference vectors to reference tensors (2D)
        PylithScalar refStressTensor[_dim * _dim];
        PylithScalar refStrainTensor[_dim * _dim];

        refStressTensor[0] = refStressVector[0];
        refStressTensor[1] = refStressVector[3];
        refStressTensor[2] = refStressVector[3];
        refStressTensor[3] = refStressVector[1];

        refStrainTensor[0] = refStrainVector[0];
        refStrainTensor[1] = refStrainVector[3];
        refStrainTensor[2] = refStrainVector[3];
        refStrainTensor[3] = refStrainVector[1];

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                g1[c * _dim + d] -= refStressTensor[c * _dim + d] + 2.0 * shearModulus * ((displacement_x[c * _dim + d] + displacement_x[d * _dim + c]) / 2.0 - refStrainTensor[c * _dim + d]);
            } // for
            g1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * (trace_strain - ref_trace_strain);
            // Biot Effective Stress Pressure Correction
            g1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // g1v_refstate

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Poroelasticity::StrainContext strainContext;
        pylith::fekernels::Poroelasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Poroelasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::PoroelasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops2D,
            stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
        const PylithInt _dim = 2;assert(_dim == dim);

        pylith::fekernels::Poroelasticity::StrainContext strainContext;
        pylith::fekernels::Poroelasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops2D);

        pylith::fekernels::Poroelasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::PoroelasticityPlaneStrain::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops2D, stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 2;
        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosityPrev = 3;

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_porosityPrev] >= 0);
        assert(porosity);

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar porosityPrev = a[aOff[i_porosityPrev]];

        // Update porosity
        porosity[0] = porosityPrev + dt * ((biotCoefficient - porosityPrev) * trace_strain_t +
                                           ((1.0 - biotCoefficient) * (biotCoefficient - porosityPrev)) /
                                           drainedBulkModulus * pressure_t);

    } // updatePorosityImplicit

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, explicit.
     */
    static inline
    void updatePorosityExplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 2;

        // Incoming solution fields.
        const PylithInt i_pressure = 1;
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosity = 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_porosity] >= 0);
        assert(porosity);

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
        std::cout << "dim:  " << dim << std::endl;
        std::cout << "numS:  " << numS << std::endl;
        std::cout << "numA:  " << numA << std::endl;
        std::cout << "sOff[0]:  " << sOff[0] << std::endl;
        std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
        std::cout << "s[0]:  " << s[0] << std::endl;
        std::cout << "aOff[0]:  " << aOff[0] << std::endl;
        std::cout << "a[0]:  " << a[0] << std::endl;
        std::cout << "t:  " << t << std::endl;
        std::cout << "x[0]:  " << x[0] << std::endl;
        std::cout << "numConstants:  " << numConstants << std::endl;
        std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

        // Do stuff
        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        // Update porosity
        porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                                  ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                                  drainedBulkModulus * pressure_t);
    } // updatePorosityExplicit

}; // IsotropicLinearPoroelasticityPlaneStrain

// ------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linear poroelasticity in 3D.
class pylith::fekernels::IsotropicLinearPoroelasticity3D {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // ================================= MMS =======================================

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and linear time MMS.
    static inline
    void f0_mms_ql_u(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(f0);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        for (PylithInt d = 0; d < _dim - 1; ++d) {
            f0[d] -= 2.0 * shearModulus - biotCoefficient * t;
        }
        f0[_dim - 1] -= 2.0 * lambda + 4.0 * shearModulus - biotCoefficient * t;
    } // f0_quadratic_linear_u

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and linear time MMS.
    static inline
    void f0_mms_ql_p(const PylithInt dim,
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

        // Incoming re-packed solution field.

        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        PylithScalar sum = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            sum += x[d];
        }
        f0[0] -= sum / biotModulus;
    } // f0_quadratic_linear_p

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and trigonometric time MMS.
    static inline
    void f0_mms_qt_u(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        for (PylithInt d = 0; d < _dim - 1; ++d) {
            f0[d] -= 2.0 * shearModulus - biotCoefficient * PetscCosReal(t);
        }
        f0[_dim - 1] -= 2.0 * lambda + 4.0 * shearModulus - biotCoefficient * PetscCosReal(t);
    } // f0_quadratic_trig_u

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and trigonometric time MMS.
    static inline
    void f0_mms_qt_p(const PylithInt dim,
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

        // Incoming re-packed solution field.

        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        PylithScalar sum = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            sum += x[d];
        }

        f0[0] += PetscSinReal(t) * sum / biotModulus;
    } // f0_quadratic_trig_p

    // ----------------------------------------------------------------------
    // f0u function for trigonometric space and linear time MMS.
    static inline
    void f0_mms_tl_u(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.

        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar lambda = drainedBulkModulus - 2.0 / 3.0 * shearModulus;

        for (PylithInt d = 0; d < _dim - 1; ++d) {
            f0[d] += PetscSqr(2. * PETSC_PI) * PetscSinReal(2. * PETSC_PI * x[d]) * (2. * shearModulus + lambda) + 2.0 * (shearModulus + lambda) - 2. * PETSC_PI * biotCoefficient * PetscSinReal(2. * PETSC_PI * x[d]) * t;
        }
        f0[_dim - 1] += PetscSqr(2. * PETSC_PI) * PetscSinReal(2. * PETSC_PI * x[_dim - 1]) * (2. * shearModulus + lambda) - 2. * PETSC_PI * biotCoefficient * PetscSinReal(2. * PETSC_PI * x[_dim - 1]) * t;
    } // f0_trig_linear_u

    // ----------------------------------------------------------------------
    // f0p function for trigonometric space and linear time MMS.
    static inline
    void f0_mms_tl_p(const PylithInt dim,
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

        // Incoming re-packed solution field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;
        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];
        const PylithScalar kappa = a[aOff[i_isotropicPermeability]] / a[aOff[i_fluidViscosity]];
        PylithScalar sum = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            sum += PetscCosReal(2. * PETSC_PI * x[d]);
        }

        f0[0] -= sum / biotModulus - 4 * PetscSqr(PETSC_PI) * kappa * sum * t;
    } // f0_quadratic_trig_p

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_explicit(const PylithInt dim,
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
        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += pressure_t / biotModulus;
    } // f0p_explicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms.
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
        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
    } // f0p_implicit

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
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
        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 4;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_body(const PylithInt dim,
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
        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 5;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_body

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 5;

        // Incoming re-packed auxiliary field.
        // Poroelasticity

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar source = a[aOff[i_source]];

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static inline
    void f0p_implicit_source_grav_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;
        const PylithInt i_source = 6;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(sOff_x[i_trace_strain] >= 0);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(aOff[i_biotModulus] >= 0);
        assert(aOff[i_source] >= 0);
        assert(f0);

        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        f0[0] += biotCoefficient * trace_strain_t;
        f0[0] += pressure_t / biotModulus;
        f0[0] -= source;
    } // f0p_implicit_source_grav_body

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
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

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(sOff);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);
        assert(sOff_x);
        assert(sOff_x[i_displacement] >= 0);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(f1);

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                f1[c * _dim + d] -= shearModulus * (displacement_x[c * _dim + d] + displacement_x[d * _dim + c]);
            } // for
            f1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * trace_strain;
            f1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // f1u

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void f1u_refstate(const PylithInt dim,
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

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];
        const PylithScalar trace_strain = s[sOff[i_trace_strain]];

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_rstress = numA - 7;
        const PylithInt i_rstrain = numA - 6;
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar *refStressVector = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy,
                                                                   // stress_yz,
        // stress_xz
        const PylithScalar *refStrainVector = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy,
                                                                   // strain_yz,
        // strain_xz

        const PylithScalar ref_trace_strain = refStrainVector[0] + refStrainVector[1] + refStrainVector[2];

        // Convert reference vectors to refrence tensors
        PylithScalar refStressTensor[_dim * _dim];
        PylithScalar refStrainTensor[_dim * _dim];
        PylithInt refTensorPos[9] = {0, 3, 5,
                                     3, 1, 4,
                                     5, 4, 2};

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                refStressTensor[i * _dim + j] = refStressVector[refTensorPos[i * _dim + j]];
                refStrainTensor[i * _dim + j] = refStrainVector[refTensorPos[i * _dim + j]];
            } // for
        } // for

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                f1[c * _dim + d] -= refStressTensor[c * _dim + d] + 2.0 * shearModulus * ((displacement_x[c * _dim + d] + displacement_x[d * _dim + c]) / 2.0 - refStrainTensor[c * _dim + d]);
            } // for
            f1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * (trace_strain - ref_trace_strain);
            // Biot Effective Stress Pressure Correction
            f1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // f1u_refstate

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d]);
        } // for
    } // f1p

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j]);
            } // for
        } // for
    } // f1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_body_force = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - body_force[d]);
        } // for

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_body_force = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j]);
            } // for
        } // for

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(aOff[i_gravityField] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
        } // for

    } // f1p_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_gravity_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_gravityField] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
            } // for
        } // for

    } // f1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void f1p_body_gravity(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluid_density = 1;
        const PylithInt i_fluid_viscosity = 2;

        const PylithInt i_body_force = 4;
        const PylithInt i_gravity_field = 5;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropic_permeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluid_density] >= 0);
        assert(aOff[i_fluid_viscosity] >= 0);
        assert(aOff[i_isotropic_permeability] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_gravity_field] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluid_density = a[aOff[i_fluid_density]];
        const PylithScalar fluid_viscosity = a[aOff[i_fluid_viscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

        const PylithScalar isotropic_permeability = a[aOff[i_isotropic_permeability]];

        for (PylithInt d = 0; d < _dim; ++d) {
            f1[d] += (isotropic_permeability / fluid_viscosity) * (pressure_x[d] - body_force[d] - fluid_density * gravity_field[d]);
        } // for

    } // f1p_body_gravity

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;

        const PylithInt i_body_force = 4;
        const PylithInt i_gravity_field = 5;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 4);
        assert(sOff_x);
        assert(sOff_x[i_pressure] >= 0);
        assert(aOff);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_body_force] >= 0);
        assert(aOff[i_gravity_field] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(f1);

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        const PylithScalar *body_force = &a[aOff[i_body_force]];
        const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                f1[i] += (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - body_force[j] - fluidDensity * gravity_field[j]);
            } // for
        } // for

    } // f1p_body_gravity_tensor_permeability

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 3;

        // Incoming solution field.

        // Incoming auxiliary fields.

        // Isotropic Linear Poroelasticity
        const PylithInt i_shearModulus = numA - 5;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 1);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_shearModulus] >= 0);
        assert(Jf3);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                Jf3[((i * _dim + i) * _dim + j) * _dim + j] -= shearModulus;
                Jf3[((i * _dim + j) * _dim + j) * _dim + i] -= shearModulus;
            }
        }

    } // Jf3uu

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;

        // Isotropic Linear Poroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf2);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] += biotCoefficient;
        } // for
    } // Jf2up

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        const PylithInt _dim = 3;

        // Isotropic Linear Poroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 5);
        assert(aOff);
        assert(aOff[i_drainedBulkModulus] >= 0);
        assert(aOff[i_shearModulus] >= 0);
        assert(Jf2);

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf2[d * _dim + d] -= drainedBulkModulus - (2.0 * shearModulus) / 3.0;
        } // for
    } // Jf2ue

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf3[]) {
        const PylithInt _dim = 3;

        // index of Incoming auxiliary fields.
        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_isotropicPermeability] >= 0);
        assert(Jf3);

        const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        for (PylithInt d = 0; d < _dim; ++d) {
            Jf3[d * _dim + d] += isotropicPermeablity / fluidViscosity;
        }
    } // Jf3pp

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                   PylithScalar Jf3[]) {
        const PylithInt _dim = 3;

        // index of Incoming auxiliary fields.
        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // Isotropic Linear Poroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_fluidViscosity] >= 0);
        assert(aOff[i_tensorPermeability] >= 0);
        assert(Jf3);

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                Jf3[i * _dim + j] += tensorPermeability[i * _dim + j] / fluidViscosity;
            } // for
        } // for
    } // Jf3pp_tensorPermeability

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(Jf0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        Jf0[0] = utshift / biotModulus;
    } // Jf0pp

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf0[]) {
        const PylithInt _dim = 3;

        // Incoming re-packed auxiliary field.
        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 2);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf0);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        Jf0[0] += utshift * biotCoefficient;
    } // Jf0pe

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0ppdot(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotModulus = numA - 2;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotModulus] >= 0);
        assert(Jf0);

        const PylithScalar biotModulus = a[aOff[i_biotModulus]];

        Jf0[0] += 1.0 / biotModulus;
    } // Jf0ppdot

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void Jf0pedot(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming re-packed auxiliary field.
        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_biotCoefficient] >= 0);
        assert(Jf0);

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        Jf0[0] += biotCoefficient;
    } // Jf0pedot

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_implicit

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
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
        const PylithInt _dim = 3;

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 4;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        // Poroelasticity
        const PylithInt i_source = 5;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_body

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 5;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static inline
    void g0p_source_grav_body(const PylithInt dim,
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

        // Incoming re-packed solution field.
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.
        const PylithInt i_source = 6;
        const PylithScalar source = a[aOff[i_source]];

        // IsotropicLinearPoroelasticity
        const PylithInt i_biotCoefficient = numA - 3;
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        g0[0] += source;
        g0[0] -= biotCoefficient * trace_strain_t;
    } // g0p_source_grav_body

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        for (PylithInt d = 0; d < _dim; ++d) {
            g1[d] -= (isotropicPermeability / fluidViscosity) * (pressure_x[d] - fluidDensity * gravityField[d]);
        } // for

    } // g1p_gravity

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_gravity_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_gravityField = 4;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];
        const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];
        const PylithScalar *gravityField = &a[aOff[i_gravityField]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                g1[i] -= (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j] - fluidDensity * gravityField[j]);
            } // for
        } // for

    } // g1p_gravity_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_isotropicPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar isotropicPermeability = a[aOff[i_isotropicPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        for (PylithInt d = 0; d < _dim; ++d) {
            g1[d] -= (isotropicPermeability / fluidViscosity) * pressure_x[d];
        } // for
    } // g1p

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static inline
    void g1p_tensor_permeability(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution field.
        const PylithInt i_pressure = 1;

        // Incoming auxiliary field.

        // Poroelasticity
        const PylithInt i_fluidViscosity = 2;

        // IsotropicLinearPoroelasticity
        const PylithInt i_tensorPermeability = numA - 1;

        const PylithScalar *pressure_x = &s_x[sOff_x[i_pressure]];

        const PylithScalar *vectorPermeability = &a[aOff[i_tensorPermeability]];
        const PylithScalar fluidViscosity = a[aOff[i_fluidViscosity]];

        PylithScalar tensorPermeability[9];
        tensorPermeability[0] = vectorPermeability[0];
        tensorPermeability[1] = vectorPermeability[3];
        tensorPermeability[2] = vectorPermeability[5];
        tensorPermeability[3] = vectorPermeability[3];
        tensorPermeability[4] = vectorPermeability[1];
        tensorPermeability[5] = vectorPermeability[4];
        tensorPermeability[6] = vectorPermeability[5];
        tensorPermeability[7] = vectorPermeability[4];
        tensorPermeability[8] = vectorPermeability[2];

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; j++) {
                g1[i] -= (tensorPermeability[i * _dim + j] / fluidViscosity) * (pressure_x[j]);
            } // for
        } // for
    } // g1p_tensor_permeability

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static inline
    void g1v(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain += displacement_x[d * _dim + d];
        }

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                g1[c * dim + d] -= shearModulus * (displacement_x[c * _dim + d] + displacement_x[d * _dim + c]);
            } // for
            g1[c * dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * trace_strain;
            g1[c * dim + c] += biotCoefficient * pressure;
        } // for
    } // g1v

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static inline
    void g1v_refstate(const PylithInt dim,
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
        const PylithInt _dim = 3;

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;

        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar pressure = s[sOff[i_pressure]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain += displacement_x[d * _dim + d];
        }

        // Incoming auxiliary fields.

        // IsotropicLinearPoroelasticity
        const PylithInt i_rstress = numA - 7;
        const PylithInt i_rstrain = numA - 6;
        const PylithInt i_shearModulus = numA - 5;
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        const PylithScalar shearModulus = a[aOff[i_shearModulus]];
        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
        const PylithScalar *refStressVector = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
        const PylithScalar *refStrainVector = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

        const PylithScalar ref_trace_strain = refStrainVector[0] + refStrainVector[1] + refStrainVector[2];
        // const PylithScalar mean_ref_stress = (refStressVector[0] + refStressVector[1] + refStressVector[2]) / 3.0;

        // Convert reference vectors to refrence tensors
        PylithScalar refStressTensor[_dim * _dim];
        PylithScalar refStrainTensor[_dim * _dim];
        PylithInt refTensorPos[9] = {0, 3, 5,
                                     3, 1, 4,
                                     5, 4, 2};

        for (PylithInt i = 0; i < _dim; ++i) {
            for (PylithInt j = 0; j < _dim; ++j) {
                refStressTensor[i * _dim + j] = refStressVector[refTensorPos[i * _dim + j]];
                refStrainTensor[i * _dim + j] = refStrainVector[refTensorPos[i * _dim + j]];
            } // for
        } // for

        for (PylithInt c = 0; c < _dim; ++c) {
            for (PylithInt d = 0; d < _dim; ++d) {
                g1[c * _dim + d] -= refStressTensor[c * _dim + d] + 2.0 * shearModulus * ((displacement_x[c * _dim + d] + displacement_x[d * _dim + c]) / 2.0 - refStrainTensor[c * _dim + d]);
            } // for
            g1[c * _dim + c] -= (drainedBulkModulus - (2.0 * shearModulus) / 3.0) * (trace_strain - ref_trace_strain);
            // Biot Effective Stress Pressure Correction
            g1[c * _dim + c] += biotCoefficient * pressure;
        } // for
    } // g1v_refstate

    // ===========================================================================================
    // Kernels for output
    // ===========================================================================================

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITHOUT reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Poroelasticity::StrainContext strainContext;
        pylith::fekernels::Poroelasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Poroelasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Poroelasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress,
            pylith::fekernels::Tensor::ops3D,
            stressVector);

    } // cauchyStress_infinitesimalStrain_asVector

    // --------------------------------------------------------------------------------------------
    /** Entry function for calculating Cauchy stress for 3D isotropic linear poroelasticity with
     * infinitesimal strain WITH a reference stress and strain.
     *
     * Used to output of Cauchy stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., biot_coefficient(1), shear_modulus(1), drained_bulk_modulus(1)]
     */
    static inline
    void cauchyStress_infinitesimalStrain_refState_asVector(const PylithInt dim,
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
        const PylithInt _dim = 3;assert(_dim == dim);

        pylith::fekernels::Poroelasticity::StrainContext strainContext;
        pylith::fekernels::Poroelasticity::setStrainContext(&strainContext, _dim, numS, sOff, sOff_x, s, s_t, s_x, x);

        pylith::fekernels::IsotropicLinearPoroelasticity::Context rheologyContext;
        pylith::fekernels::IsotropicLinearPoroelasticity::setContext_refState(
            &rheologyContext, _dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::Tensor::ops3D);

        pylith::fekernels::Poroelasticity::stress_asVector(
            strainContext, &rheologyContext,
            pylith::fekernels::Poroelasticity3D::infinitesimalStrain,
            pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refState,
            pylith::fekernels::Tensor::ops3D, stressVector);
    } // cauchyStress_infinitesimalStrain_refState_asVector

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, implicit.
     */
    static inline
    void updatePorosityImplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 3;

        // Incoming solution fields.
        const PylithInt i_pressure_t = 4;
        const PylithInt i_trace_strain_t = 5;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosity = 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_porosity] >= 0);
        assert(porosity);

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
        std::cout << "dim:  " << dim << std::endl;
        std::cout << "numS:  " << numS << std::endl;
        std::cout << "numA:  " << numA << std::endl;
        std::cout << "sOff[0]:  " << sOff[0] << std::endl;
        std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
        std::cout << "s[0]:  " << s[0] << std::endl;
        std::cout << "aOff[0]:  " << aOff[0] << std::endl;
        std::cout << "a[0]:  " << a[0] << std::endl;
        std::cout << "t:  " << t << std::endl;
        std::cout << "x[0]:  " << x[0] << std::endl;
        std::cout << "numConstants:  " << numConstants << std::endl;
        std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

        // Do stuff
        const PylithScalar pressure_t = s ? s[sOff[i_pressure_t]] : 0.0;
        const PylithScalar trace_strain_t = s ? s[sOff[i_trace_strain_t]] : 0.0;

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        // Update porosity
        porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                                  ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                                  drainedBulkModulus * pressure_t);
    } // updatePorosityImplicit

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material, explicit.
     */
    static inline
    void updatePorosityExplicit(const PylithInt dim,
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
                                PylithScalar porosity[]) {
        const PylithInt _dim = 3;

        // Incoming solution fields.
        const PylithInt i_pressure = 1;
        const PylithInt i_velocity = 2;

        // Incoming re-packed auxiliary field.

        // Poroelasticity
        const PylithInt i_porosity = 3;

        // Run Checks
        assert(_dim == dim);
        assert(numS >= 3);
        assert(numA >= 3);
        assert(aOff);
        assert(aOff[i_porosity] >= 0);
        assert(porosity);

        // IsotropicLinearPoroelasticity
        const PylithInt i_drainedBulkModulus = numA - 4;
        const PylithInt i_biotCoefficient = numA - 3;

        // Constants
        const PylithScalar dt = constants[0];

#if 0 // :DEBUG:
        std::cout << "dim:  " << dim << std::endl;
        std::cout << "numS:  " << numS << std::endl;
        std::cout << "numA:  " << numA << std::endl;
        std::cout << "sOff[0]:  " << sOff[0] << std::endl;
        std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
        std::cout << "s[0]:  " << s[0] << std::endl;
        std::cout << "aOff[0]:  " << aOff[0] << std::endl;
        std::cout << "a[0]:  " << a[0] << std::endl;
        std::cout << "t:  " << t << std::endl;
        std::cout << "x[0]:  " << x[0] << std::endl;
        std::cout << "numConstants:  " << numConstants << std::endl;
        std::cout << "porosity[0]:  " << totalStrain[0] << std::endl;
#endif

        // Do stuff
        const PylithScalar pressure_t = s_t[sOff[i_pressure]];
        const PylithScalar *velocity_x = &s_x[sOff[i_velocity]];

        const PylithScalar drainedBulkModulus = a[aOff[i_drainedBulkModulus]];
        const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

        PylithScalar trace_strain_t = 0.0;

        for (PylithInt d = 0; d < _dim; ++d) {
            trace_strain_t += velocity_x[d * _dim + d];
        }

        // Update porosity
        porosity[0] = a[aOff[i_porosity]] + dt * ((biotCoefficient - a[aOff[i_porosity]]) * trace_strain_t +
                                                  ((1.0 - biotCoefficient) * (biotCoefficient - a[aOff[i_porosity]])) /
                                                  drainedBulkModulus * pressure_t);
    } // updatePorosityExplicit

}; // IsotropicLinearPoroelasticity3D

#endif // pylith_fekernels_isotropiclinearporoelasticity_hh

// End of file
