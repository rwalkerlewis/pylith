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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformapp.

#include "ElasticityImplicitLgDeformGravData3DLinear.hh"

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_spaceDim = 3;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_cellDim = 3;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_numVertices = 4;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_numBasis = 4;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_numQuadPts = 1;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_matLabel = "elastic isotropic 3-D";

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_vertices[] = {
 -5.00000000e-01, -1.00000000e+00, -5.00000000e-01,
  2.00000000e+00, -5.00000000e-01, -4.00000000e-01,
  1.00000000e+00, -1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
};

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_cells[] = {
0,1,2,3,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_quadPts[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_quadWts[] = {
  1.33333333e+00,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_basis[] = {
  2.50000000e-01,  2.50000000e-01,  2.50000000e-01,
  2.50000000e-01,};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_fieldTIncr[] = {
  3.00000000e-01,  2.00000000e-01, -5.00000000e-01,
 -3.00000000e-01, -4.00000000e-01, -6.00000000e-01,
  2.00000000e-01,  6.00000000e-01,  3.00000000e-01,
 -6.00000000e-01, -1.00000000e-01, -3.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_fieldT[] = {
  8.00000000e-01,  1.00000000e-01, -6.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -5.00000000e-01,
  1.00000000e-01,  7.00000000e-01,  2.00000000e-01,
 -5.00000000e-01, -0.00000000e+00, -2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_fieldTmdt[] = {
  1.00000000e-01,  1.00000000e-01, -3.00000000e-01,
 -2.00000000e-01, -1.00000000e-01, -5.00000000e-01,
  2.00000000e-01,  4.00000000e-01,  1.00000000e-01,
 -4.00000000e-01, -1.00000000e-01, -1.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_valsResidual[] = {
  2.79474210e+11,  1.80174876e+12,  1.23194350e+12,
  6.44632598e+11,  4.29759382e+12,  2.85238668e+12,
 -1.13174893e+12, -6.94850554e+12, -4.68862818e+12,
  2.07642124e+11,  8.49162957e+11,  4.69797995e+11,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::_valsJacobian[] = {
  2.94874367e+11,  6.16951570e+10,  3.89461610e+10,
  5.48150595e+11,  1.11804684e+11,  7.10898629e+10,
 -9.45152856e+11, -2.08621385e+11, -1.32482083e+11,
  1.02127894e+11,  3.51215443e+10,  2.24460592e+10,
  6.16951570e+10,  6.57018660e+11,  2.71425583e+11,
  1.11804684e+11,  1.23922444e+12,  4.97844328e+11,
 -2.08621385e+11, -2.12462667e+12, -8.51891691e+11,
  3.51215443e+10,  2.28383572e+11,  8.26217801e+10,
  3.89461610e+10,  2.71425583e+11,  4.77828353e+11,
  7.10898629e+10,  4.97844328e+11,  8.73012939e+11,
 -1.32482083e+11, -8.51891691e+11, -1.49263907e+12,
  2.24460592e+10,  8.26217801e+10,  1.41797778e+11,
  5.48150595e+11,  1.11804684e+11,  7.10898629e+10,
  1.41837229e+12,  2.89096417e+11,  1.70605497e+11,
 -2.25396787e+12, -4.94966564e+11, -2.95400181e+11,
  2.87444987e+11,  9.40654624e+10,  5.37048209e+10,
  1.11804684e+11,  1.23922444e+12,  4.97844328e+11,
  2.89096417e+11,  3.27225597e+12,  1.26969767e+12,
 -4.94966564e+11, -5.16973245e+12, -1.99940566e+12,
  9.40654624e+10,  6.58252043e+11,  2.31863658e+11,
  7.10898629e+10,  4.97844328e+11,  8.73012939e+11,
  1.70605497e+11,  1.26969767e+12,  2.20481950e+12,
 -2.95400181e+11, -1.99940566e+12, -3.47904176e+12,
  5.37048209e+10,  2.31863658e+11,  4.01209317e+11,
 -9.45152856e+11, -2.08621385e+11, -1.32482083e+11,
 -2.25396787e+12, -4.94966564e+11, -2.95400181e+11,
  3.68810232e+12,  8.67069202e+11,  5.19506826e+11,
 -4.88981593e+11, -1.63481253e+11, -9.16245619e+10,
 -2.08621385e+11, -2.12462667e+12, -8.51891691e+11,
 -4.94966564e+11, -5.16973245e+12, -1.99940566e+12,
  8.67069202e+11,  8.39475093e+12,  3.23896409e+12,
 -1.63481253e+11, -1.10039182e+12, -3.87666745e+11,
 -1.32482083e+11, -8.51891691e+11, -1.49263907e+12,
 -2.95400181e+11, -1.99940566e+12, -3.47904176e+12,
  5.19506826e+11,  3.23896409e+12,  5.64558477e+12,
 -9.16245619e+10, -3.87666745e+11, -6.73903937e+11,
  1.02127894e+11,  3.51215443e+10,  2.24460592e+10,
  2.87444987e+11,  9.40654624e+10,  5.37048209e+10,
 -4.88981593e+11, -1.63481253e+11, -9.16245619e+10,
  9.94087115e+10,  3.42942460e+10,  1.54736818e+10,
  3.51215443e+10,  2.28383572e+11,  8.26217801e+10,
  9.40654624e+10,  6.58252043e+11,  2.31863658e+11,
 -1.63481253e+11, -1.10039182e+12, -3.87666745e+11,
  3.42942460e+10,  2.13756204e+11,  7.31813067e+10,
  2.24460592e+10,  8.26217801e+10,  1.41797778e+11,
  5.37048209e+10,  2.31863658e+11,  4.01209317e+11,
 -9.16245619e+10, -3.87666745e+11, -6.73903937e+11,
  1.54736818e+10,  7.31813067e+10,  1.30896842e+11,
};

pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::ElasticityImplicitLgDeformGravData3DLinear(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityImplicitLgDeformGravData3DLinear::~ElasticityImplicitLgDeformGravData3DLinear(void)
{}


// End of file
