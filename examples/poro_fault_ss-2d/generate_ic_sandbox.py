#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import numpy 
from netCDF4 import Dataset
import scipy.ndimage
import matplotlib.pyplot as plt

# mesh = Dataset('mesh_hex.exo','r')

# Dimensions
DOMAIN_X = 150.0e+3
DOMAIN_Y = 100.0e+3

# Discretization
dx = 5.0e+3
dy = 5.0e+3 

nx = numpy.int64(DOMAIN_X / dx)
ny = numpy.int64(DOMAIN_Y / dy)

# Grid
x = numpy.linspace(-DOMAIN_X/2,DOMAIN_X/2,nx+1)
y = numpy.linspace(-DOMAIN_Y/2,DOMAIN_Y/2,ny+1)

xx, yy = numpy.meshgrid(x,y)

# Initial time
t0 = 0.0


# Two Dimensional Values
u_x = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy - 1)
u_y = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy)
p = t0*(yy*yy - 1)
trace_strain = numpy.nan_to_num(yy / numpy.abs(yy), nan=0.0) * 2 * yy
L_x = numpy.ones(x.size) * 0
L_y = numpy.ones(y.size) * t0
p_f = numpy.ones(x.size) * -4.0 * 0.0 + t0

# Solution


# Aux
fluid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
porosity = 0.5 * numpy.ones(nx*ny) # kg / m**3
biot_coefficient = 1.0 * numpy.ones(nx*ny)
fluid_viscosity = 1.0 * numpy.ones(nx*ny) # Pa*s
shear_modulus = 0.5 * numpy.ones(nx*ny) # Pa
fluid_bulk_modulus = 2e9 * numpy.ones(nx*ny) # Pa
drained_bulk_modulus = 10e9 * numpy.ones(nx*ny) # Pa
solid_bulk_modulus = 11039657020.4 * numpy.ones(nx*ny) # Pa

# Convert data to array for DB generation



# End of file

