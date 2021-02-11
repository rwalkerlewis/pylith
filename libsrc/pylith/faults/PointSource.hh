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
    ~PointSource(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set point sources.
     *
     * @param names Array of point source names.
     * @param numNames Number of point source names.
     * @param ruptures Array of point sources.
     * @param numRuptures Number of point sources.
     */
    void setPointSources(const char* const* names,
                         const int numNames,
                         KinSrc** ruptures,
                         const int numRuptures);




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

    /** Set time delay of point source.
     *
     * @param vec Reference direction unit vector.
     */
    void setPointShift(const PylithReal vec[1]);

    /** Set dominant frequency of Ricker function.
     *
     * @param vec Reference direction unit vector.
     */
    void setDominantFrequency(const PylithReal vec[1]);

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

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;
    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:
    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE TYPEDEFS ////////////////////////////////////////////////////////////////////////////////////////////////
private:

    typedef std::map<std::string, PtSrc*> srcs_type;    

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _momentTensor[9]; ///< Full representation of the moment tensor.
    PylithReal _pointLocation[3]; ///< Cartesian representation of the origin of the point source.
    PylithReal _originTime[1]; ///< Elapsed time to pass for implementation of the point source.
    PylithReal _dominantFrequency[1]; ///< Dominant frequency (Hz) for Ricker source function.

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    pylith::faults::AuxiliaryFactoryKinematic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    srcs_type _ruptures; ///< Array of kinematic earthquake ruptures.
    std::string _sourceLabel; ///< Label identifying point source.
    std::string _momentTensorConvention; ///< Label identifying convention of moment tensor reference frame.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PointSource(const PointSource&); ///< Not implemented
    const PointSource& operator=(const PointSource&); ///< Not implemented

}; // class PointSource

#endif // pylith_faults_pointsource_hh

// End of file
