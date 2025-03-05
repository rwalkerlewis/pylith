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

// PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

static const int GET_SLIP;
static const int GET_SLIP_RATE;
static const int GET_SLIP_ACC;


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


    /** Get requested slip subfields at time t.
     *
     * @param[inout] slipLocalVec Local PETSc vector for slip, slip rate, or slip accelerationvalues.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] bitSlipSubfields Slip subfields to compute.
     */
    virtual
    void getSlipSubfields(PetscVec slipLocalVec,
                          pylith::topology::Field* faultAuxiliaryField,
                          const PylithScalar t,
                          const PylithScalar timeScale,
                          const int bitSlipSubfields);

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
     * - numA - 3: body_force(dim)
     * - numA - 2: source(1)
     * - numA - 1: slip(dim)
     */
    // Add other kernels
    PetscPointFunc _thicknessFnKernel; ///< Kernel for fault thickness function.
    PetscPointFunc _porosityFnKernel; ///< Kernel for porosity function.
    PetscPointFunc _beta_pFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _beta_sigmaFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _fault_permeabilityFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _fluid_viscosityFnKernel; ///< Kernel for slip time function.

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
