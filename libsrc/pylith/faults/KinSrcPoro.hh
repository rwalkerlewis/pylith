// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/KinSrcPoro.hh
 *
 * @brief C++ object for managing parameters for a kinematic
 * earthquake source.
 */

#if !defined(pylith_faults_kinsrcporo_hh)
#define pylith_faults_kinsrcporo_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::SpaceDim
#include "pylith/feassemble/feassemblefwd.hh" // USES AuxiliaryFactory
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// KinSrcPoro -------------------------------------------------------------
/** @brief Kinematic earthquake source.
 *
 * KinSrcPoro is responsible for providing the value of slip at time t
 * over a fault surface.
 *
 * The fault integrator's auxiliary field has the 'slip' subfield. The
 * auxiliary subfields in this object use the discretization from the
 * 'slip' subfield.
 */
class pylith::faults::KinSrcPoro : public pylith::utils::PyreComponent {
    friend class TestKinSrcPoro; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcPoro(void);

    /// Destructor.
    virtual ~KinSrcPoro(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set origin time for earthquake source.
     *
     * @param value Origin time for earthquake source.
     */
    void originTime(const PylithReal value);

    /** Get origin time for earthquake source.
     *
     * @returns Origin time for earthquake source.
     */
    PylithReal originTime(void) const;

    /** Get auxiliary field associated with the kinematic source.
     *
     * @return field Auxiliary field for the kinematic source.
     */
    const pylith::topology::Field& auxField(void) const;

    /** Set the spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Initialize kinematic (prescribed slip) earthquake source.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void initialize(const pylith::topology::Field& auxField,
                    const spatialdata::units::Nondimensional& normalizer,
                    const spatialdata::geocoords::CoordSys* cs);

    // ** TO DO **
    // implement the following add* functions
    /** -- numA : number of auxiliary fields
     ***** Required fields
     * - 0: thickness(1)
     * - 1: porosity(1)
     * - 2: beta_p(1)
     * - 3: beta_sigma(1)
     * - 4: permeability_tangential(1)
     * - 5: permeability_normal(1)
     * - 6: fluid_viscosity(1)
     * - 7: bulk_modulus_negative(1)
     * - 8: shear_modulus_negative(1)
     * - 9: bulk_modulus_positive(1)
     * - 10: shear_modulus_positive(1)
     * - numA - 3: body_force(dim)
     * - numA - 2: source(1)
     * - numA - 1: slip(dim)
     */

    /** Set thickness values at time t.
     *
     * @param[inout] thicknessLocalScalar Local PETSc scalar for thickness value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatethickness(PetscVec thicknessLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set porosity at time t.
     *
     * @param[inout] porosityLocalScalar Local PETSc scalar for porosity value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateporosity(PetscVec porosityLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set porosity at time t.
     *
     * @param[inout] beta_pLocalScalar Local PETSc scalar for beta_p value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatebeta_p(PetscVec beta_pLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set beta_sigma at time t.
     *
     * @param[inout] beta_sigmaLocalScalar Local PETSc scalar for beta_sigma value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatebeta_sigma(PetscVec beta_sigmaLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set permeability_tangential at time t.
     *
     * @param[inout] permeability_tangentialLocalScalar Local PETSc scalar for permeability_tangential value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatepermeability_tangential(PetscVec permeability_tangentialLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set permeability_normal at time t.
     *
     * @param[inout] permeability_normalLocalScalar Local PETSc scalar for permeability_normal value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatepermeability_normal(PetscVec permeability_normalLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set fluid_viscosity at time t.
     *
     * @param[inout] fluid_viscosityLocalScalar Local PETSc scalar for fluid_viscosity value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatefluid_viscosity(PetscVec fluid_viscosityLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set bulk_modulus_negative at time t.
     *
     * @param[inout] bulk_modulus_negativeLocalScalar Local PETSc scalar for bulk_modulus_negative value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatebulk_modulus_negative(PetscVec bulk_modulus_negativeLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set shear_modulus_negative at time t.
     *
     * @param[inout] shear_modulus_negativeLocalScalar Local PETSc scalar for shear_modulus_negative value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateshear_modulus_negative(PetscVec shear_modulus_negativeLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set bulk_modulus_positive at time t.
     *
     * @param[inout] bulk_modulus_positiveLocalScalar Local PETSc scalar for bulk_modulus_positive value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updatebulk_modulus_positive(PetscVec bulk_modulus_positiveLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);

    /** Set shear_modulus_positive at time t.
     *
     * @param[inout] shear_modulus_positiveLocalScalar Local PETSc scalar for sheaer_modulus_positive value.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateshear_modulus_positive(PetscVec shear_modulus_positiveLocalVec,
                         pylith::topology::Field* faultAuxiliaryField,
                         const PylithScalar t,
                         const PylithScalar timeScale);
    
    /** Set slip values at time t.
     *
     * @param[inout] slipLocalVec Local PETSc vector for slip values.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateSlip(PetscVec slipLocalVec,
                    pylith::topology::Field* faultAuxiliaryField,
                    const PylithScalar t,
                    const PylithScalar timeScale);

    /** Set slip rate values at time t.
     *
     * @param[inout] slipRateLocalVec Local PETSc vector for slip values.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateSlipRate(PetscVec slipRateLocalVec,
                        pylith::topology::Field* faultAuxiliaryField,
                        const PylithScalar t,
                        const PylithScalar timeScale);

    /** Set slip acceleration values at time t.
     *
     * @param[inout] slipAccLocalVec Local PETSc vector for slip values.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateSlipAcc(PetscVec slipAccLocalVec,
                       pylith::topology::Field* faultAuxiliaryField,
                       const PylithScalar t,
                       const PylithScalar timeScale);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    virtual
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                              const spatialdata::geocoords::CoordSys* cs) = 0;

    /** Set constants used in finite-element integrations.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     */
    void _setFEConstants(const pylith::topology::Field& auxField) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::faults::KinSrcPoroAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    // implement the following add* functions
    /** -- numA : number of auxiliary fields
     ***** Required fields
     * - 0: thickness(1)
     * - 1: porosity(1)
     * - 2: beta_p(1)
     * - 3: beta_sigma(1)
     * - 4: permeability_tangential(1)
     * - 5: permeability_normal(1)
     * - 6: fluid_viscosity(1)
     * - 7: bulk_modulus_negative(1)
     * - 8: shear_modulus_negative(1)
     * - 9: bulk_modulus_positive(1)
     * - 10: shear_modulus_positive(1)
     * - numA - 3: body_force(dim)
     * - numA - 2: source(1)
     * - numA - 1: slip(dim)
     */
    // Add other kernels
    PetscPointFunc _thicknessFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _porosityFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _beta_pFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _beta_sigmaFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _permeability_tangentialFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _permeability_normalFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _fluid_viscosityFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _bulk_modulus_negativeFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _shear_modulus_negativeFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _bulk_modulus_positiveFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _shear_modulus_positiveFnKernel; ///< Kernel for slip time function.

    // Original slip function kernels inherited from KinSrc
    PetscPointFunc _slipFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _slipRateFnKernel; ///< Kernel for slip rate time function.
    PetscPointFunc _slipAccFnKernel; ///< Kernel for slip acceleration time function.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for this integrator.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    PylithReal _originTime; ///< Origin time for earthquake source

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    KinSrcPoro(const KinSrcPoro&); ///< Not implemented
    const KinSrcPoro& operator=(const KinSrcPoro&); ///< Not implemented

}; // class KinSrcPoro

#endif // pylith_faults_KinSrcPoro_hh

// End of file
