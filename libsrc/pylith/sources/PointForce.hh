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

#include "pylith/sources/sourcesfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/feassemble/feassemblefwd.hh" // USES IntegratorDomain

class pylith::sources::PointForce : public pylith::problems::Physics {
    friend class TestPointForce; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    PointForce(void);

    /// Destructor.
    ~PointForce(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set source location.
     *
     * @param[in] location Source location coordinates.
     * @param[in] size Size of location array (spatial dimension).
     */
    void setLocation(const PylithReal* location,
                     const int size);

    /** Get source location.
     *
     * @returns Source location coordinates.
     */
    const PylithReal* getLocation(void) const;

    /** Get spatial dimension of source location.
     *
     * @returns Spatial dimension.
     */
    int getLocationSize(void) const;

    /** Set moment tensor components.
     *
     * The moment tensor is specified in Voigt notation:
     * - 2D: [Mxx, Myy, Mxy]
     * - 3D: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
     *
     * @param[in] momentTensor Moment tensor components.
     * @param[in] size Number of components.
     */
    void setMomentTensor(const PylithReal* momentTensor,
                         const int size);

    /** Get moment tensor components.
     *
     * @returns Moment tensor components.
     */
    const PylithReal* getMomentTensor(void) const;

    /** Get number of moment tensor components.
     *
     * @returns Number of moment tensor components.
     */
    int getMomentTensorSize(void) const;

    /** Set source magnitude (seismic moment).
     *
     * @param[in] value Magnitude value.
     */
    void setMagnitude(const PylithReal value);

    /** Get source magnitude.
     *
     * @returns Source magnitude.
     */
    PylithReal getMagnitude(void) const;

    /** Set origin time.
     *
     * @param[in] value Origin time.
     */
    void setOriginTime(const PylithReal value);

    /** Get origin time.
     *
     * @returns Origin time.
     */
    PylithReal getOriginTime(void) const;

    /** Set dominant frequency for source time function.
     *
     * @param[in] value Dominant frequency.
     */
    void setDominantFrequency(const PylithReal value);

    /** Get dominant frequency.
     *
     * @returns Dominant frequency.
     */
    PylithReal getDominantFrequency(void) const;

    /** Set time delay for source time function.
     *
     * @param[in] value Time delay.
     */
    void setTimeDelay(const PylithReal value);

    /** Get time delay.
     *
     * @returns Time delay.
     */
    PylithReal getTimeDelay(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in] domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in] domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    /** Update auxiliary field values to current time.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
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

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for source.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for source.
     * @param[in] solution Solution field.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::vector<PylithReal> _location; ///< Source location.
    std::vector<PylithReal> _momentTensor; ///< Moment tensor components.
    PylithReal _magnitude; ///< Source magnitude (seismic moment).
    PylithReal _originTime; ///< Origin time.
    PylithReal _dominantFrequency; ///< Dominant frequency for source time function.
    PylithReal _timeDelay; ///< Time delay for source time function.

    pylith::sources::PointForceAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    PointForce(const PointForce&); ///< Not implemented.
    const PointForce& operator=(const PointForce&); ///< Not implemented.

}; // class PointForce

// End of file
