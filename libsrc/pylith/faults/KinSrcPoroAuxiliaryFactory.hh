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

/** @file libsrc/faults/KinSrcPoroAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for slip time functions.
 */

#if !defined(pylith_faults_kinsrcporoauxiliaryfactory_hh)
#define pylith_faults_kinsrcporoauxiliaryfactory_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::KinSrcPoroAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestSlipFnAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcPoroAuxiliaryFactory(void);

    /// Destructor.
    virtual ~KinSrcPoroAuxiliaryFactory(void);

    /// Add slip initiation time (relative to origin time) subfield to auxiliary field.
    void addInitiationTime(void);

    /// Add rise time subfield to auxiliary field.
    void addRiseTime(void);

    /// Add _slip subfield to auxiliary field.
    void addSlip(void);

    // Newly added components
    /// Add thickness subfield to auxiliary field
    void addThickness(void);

    /// Add porosity subfield to auxiliary field
    void addPorosity(void);

    /// Add beta_p subfield to auxiliary field
    void addBetaP(void);

    /// Add beta_sigma subfield to auxiliary field
    void addBetaSigma(void);

    /// Add permeability_tangential subfield to auxiliary field
    void addPermeabilityTangential(void);

    /// Add permeability_normal subfield to auxiliary field
    void addPermeabilityNormal(void);

    /// Add fluid viscosity subfield to auxiliary field
    void addFluidViscosity(void);

    /// Add bulkmodulusnegative subfield to auxiliary field
    void addBulkModulusNegative(void);

    /// Add bulkmoduluspositive subfield to auxiliary field
    void addBulkModulusPositive(void);

    /// Add shearmodulusnegative subfield to auxiliary field
    void addShearModulusNegative(void);

    /// Add shearmoduluspositive subfield to auxiliary field
    void addShearModulusPositive(void);

    /// Add bodyforce subfield to auxiliary field
    void addBodyForce(void);

    /// Add source subfield to auxiliary field
    void addSource(void);

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

    KinSrcPoroAuxiliaryFactory(const KinSrcPoroAuxiliaryFactory&); ///< Not implemented.
    const KinSrcPoroAuxiliaryFactory& operator=(const KinSrcPoroAuxiliaryFactory&); ///< Not implemented

}; // class KinSrcPoroAuxiliaryFactory

#endif // pylith_faults_kinsrcporoauxiliaryfactory_hh

// End of file
