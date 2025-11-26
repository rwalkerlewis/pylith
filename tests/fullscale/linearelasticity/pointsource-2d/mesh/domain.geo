// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//
// Gmsh geometry file for 2D point source test.
// Domain: -8 km <= x <= 8 km, -8 km <= y <= 8 km
// Units: km

// Parameters
domain_x = 8.0;
domain_y = 8.0;
h = 1.0;  // element size

// Points
Point(1) = {-domain_x, -domain_y, 0, h};
Point(2) = { domain_x, -domain_y, 0, h};
Point(3) = { domain_x,  domain_y, 0, h};
Point(4) = {-domain_x,  domain_y, 0, h};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("elastic", 1) = {1};
Physical Curve("boundary_yneg", 12) = {1};
Physical Curve("boundary_xpos", 11) = {2};
Physical Curve("boundary_ypos", 13) = {3};
Physical Curve("boundary_xneg", 10) = {4};


// End of file
