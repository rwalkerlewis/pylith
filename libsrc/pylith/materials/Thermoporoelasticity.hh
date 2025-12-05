// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

class pylith::materials::Thermoporoelasticity : public pylith::materials::Material {
    friend class TestThermoporoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Thermoporoelasticity(void);

    /// Destructor.
    ~Thermoporoelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Include body force?
     *
     * @param[in] value Flag indicating to include body force term.
     */
    void useBodyForce(const bool value);

    /** Include body force?
     *
     * @returns True if including body force term, false otherwise.
     */
    bool useBodyForce(void) const;

    /** Include source density (fluid source)?
     *
     * @param[in] value Flag indicating to include source density term.
     */
    void useSourceDensity(const bool value);

    /** Include source density?
     *
     * @returns True if including source density term, false otherwise.
     */
    bool useSourceDensity(void) const;

    /** Include heat source?
     *
     * @param[in] value Flag indicating to include heat source term.
     */
    void useHeatSource(const bool value);

    /** Include heat source?
     *
     * @returns True if including heat source term, false otherwise.
     */
    bool useHeatSource(void) const;

    /** Use reference stress and strain in computation of stress and strain?
     *
     * @param[in] value Flag indicating to include reference stress and strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Set bulk rheology.
     *
     * @param[in] rheology Bulk rheology for thermoporoelasticity.
     */
    void setBulkRheology(pylith::materials::RheologyThermoporoelasticity* const rheology);

    /** Get bulk rheology.
     *
     * @returns Bulk rheology for thermoporoelasticity.
     */
    pylith::materials::RheologyThermoporoelasticity* getBulkRheology(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     *
     *  @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    /** Get default PETSc solver options appropriate for material.
     *
     * @param[in] isParallel True if running in parallel, False if running in serial.
     * @param[in] hasFault True if problem has fault, False otherwise.
     * @returns PETSc solver options.
     */
    pylith::utils::PetscOptions* getSolverDefaults(const bool isParallel,
                                                   const bool hasFault) const;

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Get derived factory associated with physics.
     *
     * @return Derived factory for physics object.
     */
    pylith::topology::FieldFactory* _getDerivedFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    /** Set kernels for computing derived field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                 const pylith::topology::Field& solution) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    bool _useBodyForce; ///< Flag to include body force term.
    bool _useSourceDensity; ///< Flag to include fluid source density term.
    bool _useHeatSource; ///< Flag to include heat source term.
    bool _useReferenceState; ///< Flag to use reference stress and strain.
    pylith::materials::RheologyThermoporoelasticity* _rheology; ///< Bulk rheology for thermoporoelasticity.
    pylith::materials::DerivedFactoryPoroelasticity* _derivedFactory; ///< Factory for creating derived fields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Thermoporoelasticity(const Thermoporoelasticity&); ///< Not implemented.
    const Thermoporoelasticity& operator=(const Thermoporoelasticity&); /// Not implemented.

}; // class Thermoporoelasticity

// End of file
