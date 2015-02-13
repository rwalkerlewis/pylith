// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//


/* Original mesh
 *
 * Cells are 0-1 and vertices are 2-13.
 *
 *       2,3,4,5 -------- 6,7,8,9 -------- 10,11,12,13
 *
 *                        ^^^^^^^ Vertices forming fault
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 and vertices are 3-18,19-22.
 *
 *       3,4,5,6 -------- 7,8,9,10 -- 15,16,17,18 -------- 11,12,13,14
 *                                    59,60,61,62
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveDynDataHex8.hh"

const char* pylith::faults::CohesiveDynDataHex8::_meshFilename =
  "data/hex8.mesh";

const int pylith::faults::CohesiveDynDataHex8::_spaceDim = 3;

const int pylith::faults::CohesiveDynDataHex8::_cellDim = 2;

const int pylith::faults::CohesiveDynDataHex8::_numBasis = 4;

const int pylith::faults::CohesiveDynDataHex8::_numQuadPts = 4;

const PylithScalar pylith::faults::CohesiveDynDataHex8::_quadPts[] = {
  -1.0, -1.0,
  +1.0, -1.0,
  +1.0, +1.0,
  -1.0, +1.0
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_quadWts[] = {
  1.0, 1.0, 1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_basis[] = {
  1.0, 0.0, 0.0, 0.0,
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_basisDeriv[] = {
  -0.39433757, -0.39433757,
  +0.39433757, -0.10566243,
  +0.10566243, +0.10566243,
  -0.10566243, +0.39433757,

  -0.39433757, -0.10566243,
  +0.39433757, -0.39433757,
  +0.10566243, +0.39433757,
  -0.10566243, +0.10566243,

  -0.10566243, -0.10566243,
  +0.10566243, -0.39433757,
  +0.39433757, +0.39433757,
  -0.39433757, +0.10566243,

  -0.10566243, -0.39433757,
  +0.10566243, -0.10566243,
  +0.39433757, +0.10566243,
  -0.39433757, +0.39433757,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_verticesRef[] = {
  -1.0, -1.0,
  +1.0, -1.0,
  +1.0, +1.0,
  -1.0, +1.0
};

const int pylith::faults::CohesiveDynDataHex8::_id = 10;

const char* pylith::faults::CohesiveDynDataHex8::_label = "fault";

const char* pylith::faults::CohesiveDynDataHex8::_initialTractFilename = 
  "data/hex8_initialtract.spatialdb";

const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldT[] = {
  4.1, 2.1, 3.1,
  4.2, 2.2, 3.2,
  4.3, 2.3, 3.3,
  4.4, 2.4, 3.4,
  4.5, 2.5, 3.5, // 6
  4.6, 2.6, 3.6, // 7
  4.7, 2.7, 3.7, // 8
  4.8, 2.8, 3.8, // 9
  4.9, 2.9, 3.9,
  4.0, 2.0, 3.0,
  4.1, 2.1, 3.1,
  4.2, 2.2, 3.2,
  4.5, 3.2, 4.3, // 15
  4.6, 3.5, 4.6, // 16
  4.7, 3.7, 4.6, // 17
  4.8, 3.6, 4.5, // 18
 -4.4, 2.4, 3.4, // 59
 -4.6, 2.6, 3.6, // 60
 -4.8, 2.8, 3.8, // 61
 -4.0, 2.0, 3.0, // 62
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_jacobian[] = {
  1.0,  0.1,  0.2, // 2x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 2y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 2z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 3x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 3y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 3z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 4x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 4y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 4z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 5x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 5y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 5z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 6x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
 +6.0, -0.5, -0.6, // 6
 -0.7, -0.8, -0.9, // 7
 -1.0, -0.8, -0.7, // 8
 -0.6, -0.5, -0.4, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 6y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.5, +6.1, -1.0, // 6
 -1.1, -1.2, -1.3, // 7
 -1.4, -1.3, -1.2, // 8
 -1.1, -1.0, -0.9, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 6z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.6, -1.0, +6.2, // 6
 -0.5, -0.6, -0.7, // 7
 -0.8, -0.9, -0.8, // 8
 -0.7, -0.6, -0.5, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 7x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
 -0.7, -1.1, -0.5, // 6
 +6.3, -0.8, -0.7, // 7
 -0.6, -0.5, -0.4, // 8
 -0.3, -0.2, -0.1, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 7y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.8, -1.2, -0.6, // 6
 -0.8, +6.4, -0.3, // 7
 -0.4, -0.5, -0.6, // 8
 -0.7, -0.8, -0.9, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 7z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.9, -1.3, -0.7, // 6
 -0.7, -0.3, +6.5, // 7
 -0.3, -0.8, -0.7, // 8
 -0.6, -0.9, -0.7, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 8x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
 -1.0, -1.4, -0.8, // 6
 -0.6, -0.4, -0.3, // 7
 +6.6, -1.1, -0.8, // 8
 -0.7, -0.6, -0.5, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 8y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.8, -1.3, -0.9, // 6
 -0.5, -0.5, -0.8, // 7
 -1.1, +6.7, -0.8, // 8
 -0.9, -1.0, -1.1, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 8z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.7, -1.2, -0.8, // 6
 -0.4, -0.6, -0.7, // 7
 -0.8, -0.8, +6.8, // 8
 -1.0, -1.1, -1.2, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 9x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
 -0.6, -1.1, -0.7, // 6
 -0.3, -0.7, -0.6, // 7
 -0.7, -0.9, -1.0, // 8
 +6.9, -0.5, -0.4, // 9
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 9y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.5, -1.0, -0.6, // 6
 -0.2, -0.8, -0.9, // 7
 -0.6, -1.0, -1.1, // 8
 -0.5, +6.0, -1.2, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 9z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
 -0.4, -0.9, -0.5, // 6
 -0.1, -0.9, -0.7, // 7
 -0.5, -1.1, -1.2, // 8
 -0.4, -1.2, +6.1, // 9
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 10x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 10y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 10z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 11x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 11y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 11z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 12x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 12y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 12z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 13x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  1.2,  1.3,  1.4,
  1.5,  1.6,  1.7,
  1.8,  1.9,  2.0,
  2.1,  2.2,  2.3,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 13y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 13z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  3.2,  5.3,  3.4,
  1.5,  4.6,  2.7,
  2.8,  3.9,  2.0,
  2.1,  6.2,  2.3,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 14x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
 +7.0, -0.5, -0.6, // 14
 -0.7, -0.8, -0.9, // 15
 -1.0, -0.8, -0.7, // 16
 -0.6, -0.5, -0.4, // 17
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 14y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.5, +7.1, -1.0, // 14
 -1.1, -1.2, -1.3, // 15
 -1.4, -1.3, -1.2, // 16
 -1.1, -1.0, -0.9, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 14z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.6, -1.0, +7.2, // 14
 -0.5, -0.6, -0.7, // 15
 -0.8, -0.9, -0.8, // 16
 -0.7, -0.6, -0.5, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 15x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
 -0.7, -1.1, -0.5, // 14
 +7.3, -0.8, -0.7, // 15
 -0.6, -0.5, -0.4, // 16
 -0.3, -0.2, -0.1, // 17
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 15y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.8, -1.2, -0.6, // 14
 -0.8, +7.4, -0.3, // 15
 -0.4, -0.5, -0.6, // 16
 -0.7, -0.8, -0.9, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 15z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.9, -1.3, -0.7, // 14
 -0.7, -0.3, +7.5, // 15
 -0.3, -0.8, -0.7, // 16
 -0.6, -0.9, -0.7, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 16x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
 -1.0, -1.4, -0.8, // 14
 -0.6, -0.4, -0.3, // 15
 +7.6, -1.1, -0.8, // 16
 -0.7, -0.6, -0.5, // 17
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 16y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.8, -1.3, -0.9, // 14
 -0.5, -0.5, -0.8, // 15
 -1.1, +7.7, -0.8, // 16
 -0.9, -1.0, -1.1, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 16z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.7, -1.2, -0.8, // 14
 -0.4, -0.6, -0.7, // 15
 -0.8, -0.8, +7.8, // 16
 -1.0, -1.1, -1.2, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  1.0,  0.1,  0.2, // 17x
  0.3,  0.4,  0.5,
  0.6,  0.7,  0.8,
  0.9,  1.0,  1.1,
  2.4,  2.5,  2.6,
  2.7,  2.8,  2.9,
  3.0,  3.1,  3.2,
  3.3,  3.4,  3.5,
  3.6,  3.7,  3.8,
  3.9,  4.0,  4.1,
  4.2,  4.3,  4.4,
  4.5,  4.6,  4.7,
 -0.6, -1.1, -0.7, // 14
 -0.3, -0.7, -0.6, // 15
 -0.7, -0.9, -1.0, // 16
 +7.9, -0.5, -0.4, // 17
  4.8,  4.9,  5.0,
  5.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  5.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 17y
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.5, -1.0, -0.6, // 14
 -0.2, -0.8, -0.9, // 15
 -0.6, -1.0, -1.1, // 16
 -0.5, +7.0, -1.2, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  3.0,  4.1,  1.2, // 17z
  3.3,  6.4,  4.5,
  3.6,  3.7,  6.8,
  3.9,  7.0,  6.1,
  1.4,  6.5,  2.6,
  1.7,  9.8,  2.9,
  1.0,  5.1,  3.2,
  1.3,  3.4,  3.5,
  2.6,  2.7,  3.8,
  2.9,  4.0,  4.1,
  2.2,  4.3,  4.4,
  3.5,  4.6,  4.7,
 -0.4, -0.9, -0.5, // 14
 -0.1, -0.9, -0.7, // 15
 -0.5, -1.1, -1.2, // 16
 -0.4, -1.2, +7.1, // 17
  3.8,  4.9,  5.0,
  4.1,  5.2,  5.3,
  5.4,  5.5,  5.6,
  7.7,  5.8,  5.9,
  0.0,  0.0,  0.0, // 18x
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 18y
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 18z
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 19x
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 19y
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 19z
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 20x
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 20y
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 20z
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 21x
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 21y
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0, // 21z
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_orientation[] = {
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_area[] = {
  1.0, 1.0, 1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_initialTractions[] = {
  // Fault coordinate frame
  +1.0, +2.0, -3.0,
  +1.1, +2.1, -3.1,
  +1.2, +2.2, -3.2,
  +1.3, +2.3, -3.3,
};


const int pylith::faults::CohesiveDynDataHex8::_numConstraintEdges = 4;
const int pylith::faults::CohesiveDynDataHex8::_constraintEdges[] = {
  59, 60, 61, 62
};
const int pylith::faults::CohesiveDynDataHex8::_negativeVertices[] = {
   7,  8,  9, 10
};
// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldIncrStick[5*4*3] = {
  0.1, 2.1, 1.1,
  0.2, 2.2, 1.2,
  0.3, 2.3, 1.3,
  0.4, 2.4, 1.4,
  0.5, 2.5, 1.5, // 7 
  0.6, 2.6, 1.6, // 8 
  0.7, 2.7, 1.7, // 9 
  0.8, 2.8, 1.8, // 10
  0.9, 2.9, 1.9,
  0.0, 2.0, 1.0,
  1.1, 3.1, 2.1,
  1.2, 3.2, 2.2,
  0.5, 2.5, 1.5, // 15
  0.6, 2.6, 1.6, // 16
  0.7, 2.7, 1.7, // 17
  0.8, 2.8, 1.8, // 18
 -12.266666666667,  3.6, 4.6, // 59
 -12.066666666667,  5.4, 2.4, // 60
 -16.866666666667,  2.2, 8.2, // 61
 -17.666666666667, 10.0, 2.0, // 62
};

// No slip
const PylithScalar pylith::faults::CohesiveDynDataHex8::_slipStickE[] = {
  0.7,  0.8,  0.0,
  0.9,  1.0,  0.0,
  1.0,  0.9,  0.0,
  0.8,  0.7,  0.0,
};

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldIncrSlip[5*4*3] = {
  1.1, 2.1, 0.1,
  1.2, 2.2, 0.2,
  1.3, 2.3, 0.3,
  1.4, 2.4, 0.4,
  1.5, 2.5, 0.5, // 7
  1.6, 2.6, 0.6, // 8
  1.7, 2.7, 0.7, // 9
  1.8, 2.8, 0.8, // 10
  1.9, 2.9, 0.9,
  1.0, 2.0, 0.0,
  1.1, 2.1, 0.1,
  1.2, 2.2, 0.2,
  1.5, 2.5, 0.5, // 15
  1.6, 2.6, 0.6, // 16
  1.7, 2.7, 0.7, // 17
  1.8, 2.8, 0.8, // 18
 -1.4, 2.4, 0.4, // 59
 -1.6, 2.6, 0.6, // 60
 -1.8, 2.8, 0.8, // 61
 -1.0, 2.0, 0.2, // 62
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldIncrSlipE[20*3] = {
   1.100000000000,   2.100000000000,   0.100000000000,
   1.200000000000,   2.200000000000,   0.200000000000,
   1.300000000000,   2.300000000000,   0.300000000000,
   1.400000000000,   2.400000000000,   0.400000000000,
   1.500000000000,   2.500001245294,   0.500000812227,
   1.600000000000,   2.600000694980,   0.600000786664,
   1.700000000000,   2.700000891763,   0.700000882544,
   1.800000000000,   2.800000903108,   0.800000871357,
   1.900000000000,   2.900000000000,   0.900000000000,
   1.000000000000,   2.000000000000,   0.000000000000,
   1.100000000000,   2.100000000000,   0.100000000000,
   1.200000000000,   2.200000000000,   0.200000000000,
   1.500000000000,   2.499998754706,   0.499999187773,
   1.600000000000,   2.599999305020,   0.599999213336,
   1.700000000000,   2.699999108237,   0.699999117456,
   1.800000000000,   2.799999096892,   0.799999128643,
  -1.400000000000,   0.328479469127,  -1.239953753608,
  -1.600000000000,   0.293941190358,  -1.262585961634,
  -1.800000000000,   0.259995967920,  -1.286431883494,
  -1.000000000000,   0.342606428329,  -1.125914857337,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_slipSlipE[] = {
   0.699997509413,   0.799998375546,   0.000000000000,
   0.899998610039,   0.999998426672,   0.000000000000,
   0.999998216473,   0.899998234913,   0.000000000000,
   0.799998193784,   0.699998257287,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldIncrOpen[] = {
  1.1, 2.1, 0.1,
  1.2, 2.2, 0.2,
  1.3, 2.3, 0.3,
  1.4, 2.4, 0.4,
  1.5, 2.5, 0.5, // 7
  1.6, 2.6, 0.6, // 8
  1.7, 2.7, 0.7, // 9
  1.8, 2.8, 0.8, // 10
  1.9, 2.9, 0.9,
  1.0, 2.0, 0.0,
  1.1, 2.1, 0.1,
  1.2, 2.2, 0.2,
  1.5, 2.5, 0.5, // 15
  1.6, 2.6, 0.6, // 16
  1.7, 2.7, 0.7, // 17
  1.8, 2.8, 0.8, // 18
 +20.4, 2.4, 0.4, // 59
 +20.6, 2.6, 0.6, // 60
 +20.8, 2.8, 0.8, // 61
 +20.0, 2.0, 0.2, // 62
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataHex8::_fieldIncrOpenE[] = {
   1.100000000000,   2.100000000000,   0.100000000000,
   1.200000000000,   2.200000000000,   0.200000000000,
   1.300000000000,   2.300000000000,   0.300000000000,
   1.400000000000,   2.400000000000,   0.400000000000,
   1.500000000000,   2.500007882076,   0.500005293315,
   1.600000000000,   2.600004874790,   0.600005290605,
   1.700000000000,   2.700006088248,   0.700005959177,
   1.800000000000,   2.800006286383,   0.800006065901,
   1.900000000000,   2.900000000000,   0.900000000000,
   1.000000000000,   2.000000000000,   0.000000000000,
   1.100000000000,   2.100000000000,   0.100000000000,
   1.200000000000,   2.200000000000,   0.200000000000,
   1.500000000000,   2.499992117924,   0.499994706685,
   1.600000000000,   2.599995125210,   0.599994709395,
   1.700000000000,   2.699993911752,   0.699994040823,
   1.800000000000,   2.799993713617,   0.799993934099,
   4.400000000000,  -2.400000000000,  -3.400000000000,
   4.600000000000,  -2.600000000000,  -3.600000000000,
   4.800000000000,  -2.800000000000,  -3.800000000000,
   4.000000000000,  -2.000000000000,  -3.000000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataHex8::_slipOpenE[] = {
   0.699984235848,   0.799989413370,   0.000000000000,
   0.899990250419,   0.999989418790,   0.000000000000,
   0.999987823504,   0.899988081646,   0.000000000000,
   0.799987427233,   0.699987868197,   0.000000000000,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataHex8::CohesiveDynDataHex8(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  initialTractFilename = const_cast<char*>(_initialTractFilename);

  fieldT = const_cast<PylithScalar*>(_fieldT);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  initialTractions = const_cast<PylithScalar*>(_initialTractions);

  constraintEdges = const_cast<int*>(_constraintEdges);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintEdges = _numConstraintEdges;  

  // Stick
  fieldIncrStick = const_cast<PylithScalar*>(_fieldIncrStick);
  slipStickE = const_cast<PylithScalar*>(_slipStickE);

  // Slip
  fieldIncrSlip = const_cast<PylithScalar*>(_fieldIncrSlip);
  fieldIncrSlipE = const_cast<PylithScalar*>(_fieldIncrSlipE);
  slipSlipE = const_cast<PylithScalar*>(_slipSlipE);

  // Open
  fieldIncrOpen = const_cast<PylithScalar*>(_fieldIncrOpen);
  fieldIncrOpenE = const_cast<PylithScalar*>(_fieldIncrOpenE);
  slipOpenE = const_cast<PylithScalar*>(_slipOpenE);
} // constructor

pylith::faults::CohesiveDynDataHex8::~CohesiveDynDataHex8(void)
{}


// End of file
