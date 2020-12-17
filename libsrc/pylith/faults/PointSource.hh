// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/PointSource.hh
 *
 * @brief C++ abstract base class for a moment tensor point source
 *
 */

#if !defined(pylith_faults_pointsource_hh)
#define pylith_faults_pointsource_hh

#include "faultsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

class pylith::faults::PointSource : public pylith::problems::Physics {
    friend class TestPointSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    PointSource(void);

    /// Destructor.
    virtual ~PointSource(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set user identified full moment tensor.
     *
     * @param vec Reference direction unit vector.
     */
    void setMomentTensor(const PylithReal vec[9]);

    /** Set origin for point source.
     *
     * @param vec Reference direction unit vector.
     */
    void setPointLocation(const PylithReal vec[3]);

    /** Adjust mesh topology for point source implementation.
     *
     * @param mesh[in] PETSc mesh.
     */
    void adjustTopology(pylith::topology::Mesh* const mesh);

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
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

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

    /** Update auxiliary subfields at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _momentTensor[9]; ///< Full representation of the moment tensor.
    PylithReal _pointLocation[3]; ///< Cartesian representation of the origin of the point source.
    PylithReal _pointElapsed[1]; ///< Elapsed time to pass for implementation of the point source.

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int _interfaceId; ///< Identifier for cohesive cells.
    std::string _interfaceLabel; ///< Label identifying vertices associated with fault.
    std::string _buriedEdgesLabel; ///< Label identifying vertices along buried edges of fault.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PointSource(const PointSource&); ///< Not implemented
    const PointSource& operator=(const PointSource&); ///< Not implemented

}; // class PointSource

#endif // pylith_faults_pointsource_hh

// End of file
