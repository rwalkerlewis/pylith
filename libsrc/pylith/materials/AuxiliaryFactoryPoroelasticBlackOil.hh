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

/** @file libsrc/materials/AuxiliaryFactoryPoroelasticBlackOil.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for poroelastic materials.
 */

#if !defined(pylith_materials_auxiliaryfactoryporoelasticblackoil_hh)
#define pylith_materials_auxiliaryfactoryporoelasticblackoil_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/materials/AuxiliaryFactoryPoroelasticBlackOility.hh" // ISA AuxiliaryFactoryPoroelasticBlackOility

class pylith::materials::AuxiliaryFactoryPoroelasticBlackOil : public pylith::materials::AuxiliaryFactoryPoroelasticBlackOility {
    friend class TestAuxiliaryFactoryPoroelasticBlackOil; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryPoroelasticBlackOil(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryPoroelasticBlackOil(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add tensor permeability subfield to auxiliary subfields.
    void addTensorPermeability(void);

    /// Add solid Bulk Modulus subfield to auxiliary subfields.
    void addSolidBulkModulus(void);

    /// Add fluid Biot Coefficient subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add drained Bulk Modulus subfield to auxiliary subfields.
    void addDrainedBulkModulus(void);

    /// Add undrained Bulk Modulus subfield to auxiliary subfields.
    void addUndrainedBulkModulus(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add three phase fluid Bulk Modulus subfield to auxiliary subfields.
    void addThreePhaseFluidModulus(void);

    /// Add three phase fluid saturation subfield to auxiliary subfields.
    void addThreePhaseSaturation(void);

    /// Add three phase fluid viscosity subfield to auxiliary subfields.
    void addThreePhaseFluidViscosity(void);

    /// Add relative permeability subfield to auxiliary subfields.
    void addRelativePermeability(void);

    /// Add formation volume factor subfield to auxiliary subfields.
    void addFormationVolumeFactors(void);

    /// Add fluid Biot Modulus subfield to auxiliary subfields.
    void addBiotModulus(void);

    /// Add solution gas oil ratio subfield to auxiliary subfields.
    void addSolutionGasOilRatio(void);

    /// Add solution oil gas ratio subfield to auxiliary subfields.
    void addSolutionOilGasRatio(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensities(void);

    /// Add reference stress subfield to auxiliary fields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary fields.
    void addReferenceStrain(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addYoungsModulus(void);

    /// Add solid bulk modulus subfield to auxiliary subfields.
    void addPoissonsRatio(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelasticBlackOil(const AuxiliaryFactoryPoroelasticBlackOil &); ///< Not implemented.
    const AuxiliaryFactoryPoroelasticBlackOil& operator=(const AuxiliaryFactoryPoroelasticBlackOil&); ///< Not

    ///< implemented

}; // class AuxiliaryFactoryPoroelasticBlackOil

#endif // pylith_materials_auxiliaryfactoryporoelasticblackoil_hh

// End of file
