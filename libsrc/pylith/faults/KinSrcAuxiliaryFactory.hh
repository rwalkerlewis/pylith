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

/** @file libsrc/faults/KinSrcAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for slip time functions.
 */

#if !defined(pylith_faults_kinsrcauxiliaryfactory_hh)
#define pylith_faults_kinsrcauxiliaryfactory_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::KinSrcAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestSlipFnAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcAuxiliaryFactory(void);

    /// Destructor.
    virtual ~KinSrcAuxiliaryFactory(void);

    /// Add slip initiation time (relative to origin time) subfield to auxiliary field.
    void addInitiationTime(void);

    /// Add rise time subfield to auxiliary field.
    void addRiseTime(void);

    /// Add final_slip subfield to auxiliary field.
    void addFinalSlip(void);

    // Newly added components
    /// Add thickness subfield to auxiliary field
    void addFinalThickness(void);

    /// Add porosity subfield to auxiliary field
    void addFinalPorosity(void);

    /// Add beta_p subfield to auxiliary field
    void addFinalBetaP(void);

    /// Add beta_sigma subfield to auxiliary field
    void addFinalBetaSigma(void);

    /// Add permeability_tangential subfield to auxiliary field
    void addFinalPermeabilityTangential(void);

    /// Add permeability_normal subfield to auxiliary field
    void addFinalPermeabilityNormal(void);

    /// Add fluid_viscosity subfield to auxiliary field
    void addFinalFluidViscosity(void);

    /// Add bulkmodulusnegative subfield to auxiliary field
    void addFinalBulkModulusNegative(void);

    /// Add bulkmoduluspositive subfield to auxiliary field
    void addFinalBulkModulusPositive(void);

    /// Add shearmodulusnegative subfield to auxiliary field
    void addFinalShearModulusNegative(void);

    /// Add shearmoduluspositive subfield to auxiliary field
    void addFinalShearModulusPositive(void);

    /// Add bodyforce subfield to auxiliary field
    void addFinalBodyForce(void);

    /// Add source subfield to auxiliary field
    void addFinalSource(void);

    /// Add slip_rate subfield to auxiliary field.
    void addSlipRate(void);

    /// Add time_history_value subfield to auxiliary field.
    void addTimeHistoryValue(void);

    /** Update time history value subfield for current time.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current time.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] dbTimeHistory Time history database.
     */
    static
    void updateTimeHistoryValue(pylith::topology::Field* auxiliaryField,
                                const PylithReal t,
                                const PylithReal timeScale,
                                spatialdata::spatialdb::TimeHistory* const dbTimeHistory);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    KinSrcAuxiliaryFactory(const KinSrcAuxiliaryFactory&); ///< Not implemented.
    const KinSrcAuxiliaryFactory& operator=(const KinSrcAuxiliaryFactory&); ///< Not implemented

}; // class KinSrcAuxiliaryFactory

#endif // pylith_faults_kinsrcauxiliaryfactory_hh

// End of file
