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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/Source.hh
 *
 * @brief C++ abstract base class for sources.
 */

#if !defined(pylith_sources_source_hh)
#define pylith_sources_source_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

// Source -------------------------------------------------------------
/** @brief C++ abstract base class for sources.
 *
 */

class pylith::sources::Source : public pylith::problems::Physics {
    friend class TestSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Source(void);

    /// Destructor.
    virtual ~Source(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set value of label source-id used to identify source cells.
     *
     * @param value Source identifier
     */
    void setSourceId(const int value);

    /** Get value of label source-id used to identify source cells.
     *
     * @returns Source identifier
     */
    int getSourceId(void) const;

    /** Set descriptive label for source.
     *
     * @param value Label of source.
     */
    void setDescriptiveLabel(const char* value);

    /** Get descruptive label of source.
     *
     * @returns Label of source
     */
    const char* getDescriptiveLabel(void) const;

    /** Set gravity field.
     *
     * @param g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int _sourceId; ///< Value of source-id label in mesh.
    std::string _descriptiveLabel; ///< Descriptive label for source.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::scalar_array _pointCoords; ///< Array of point coordinates.
    PylithInt _numPoints; ///< Number of point coordinates
    pylith::string_vector _pointNames; ///< Array of point names.
    PylithInt _numPointsNames; ///< Number of point names

    Source(const Source&); ///< Not implemented.
    const Source& operator=(const Source&); ///< Not implemented

}; // Source

#endif // pylith_sources_source_hh

// End of file
