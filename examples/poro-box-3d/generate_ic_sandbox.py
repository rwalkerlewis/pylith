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
import matplotlib.pyplot as plt


# Dimensions
DOMAIN_X = 10.0
DOMAIN_Y = 10.0
DOMAIN_Z = 10.0

# Discretization
dx = 0.5
dy = 0.5 
dz = 0.5

nx = numpy.int64(DOMAIN_X / dx)
ny = numpy.int64(DOMAIN_Y / dy)
nz = numpy.int64(DOMAIN_Z / dz)

# Grid
x = numpy.linspace(-DOMAIN_X/2,DOMAIN_X/2,nx+1)
y = numpy.linspace(-DOMAIN_Y/2,DOMAIN_Y/2,ny+1)
z = numpy.linspace(-DOMAIN_Z/2,DOMAIN_Z/2,nz+1)

xxx, yyy, zzz = numpy.meshgrid(x,y,z)
xyz = numpy.column_stack((xxx.flatten(), yyy.flatten(), zzz.flatten()))

# Earth Tide Values

# Time in decimal days
dt = 0.01
elapsed = 3.0
tsteps = numpy.int32(elapsed/dt)
t = numpy.linspace(-dt,elapsed,tsteps+2)

def et_disp_x(x, t):
    return x*1E-8*(5*numpy.cos(t*2*numpy.pi) + 2*numpy.cos((t-0.5)*2*numpy.pi) + 1*numpy.cos((t+0.3)*0.5*numpy.pi))

def et_disp_y(y, t):
    return y*1E-8*(7*numpy.cos(t*2*numpy.pi) + 4*numpy.cos((t-0.3)*2*numpy.pi) + 7*numpy.cos((t+0.6)*0.5*numpy.pi))

def et_disp_z(z, t):
    return z*1E-8*(7*numpy.cos((t-0.5)*2*numpy.pi) + 4*numpy.cos((t-0.8)*2*numpy.pi) + 7*numpy.cos((t+0.1)*4*numpy.pi))

et_x = et_disp_x(x, t.reshape([t.size, 1]))
et_y = et_disp_y(y, t.reshape([t.size, 1]))
et_z = et_disp_z(z, t.reshape([t.size, 1]))


# # Initial time
# t0 = 0.0


# # Two Dimensional Values
# u_x = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy - 1)
# u_y = numpy.nan_to_num(yy / -numpy.abs(yy), nan=0.0) *(yy*yy)
# p = t0*(yy*yy - 1)
# trace_strain = numpy.nan_to_num(yy / numpy.abs(yy), nan=0.0) * -2 * yy
# L_x = numpy.ones(x.size) * 0
# L_y = numpy.ones(y.size) * t0
# p_f = numpy.ones(x.size) * 0.5 * t0


# # Aux
# fluid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
# solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
# solid_density = 1.0 * numpy.ones(nx*ny) # kg / m**3
# porosity = 0.5 * numpy.ones(nx*ny) # kg / m**3
# biot_coefficient = 1.0 * numpy.ones(nx*ny)
# fluid_viscosity = 1.0 * numpy.ones(nx*ny) # Pa*s
# shear_modulus = 0.5 * numpy.ones(nx*ny) # Pa
# fluid_bulk_modulus = 2e9 * numpy.ones(nx*ny) # Pa
# drained_bulk_modulus = 10e9 * numpy.ones(nx*ny) # Pa
# solid_bulk_modulus = 11039657020.4 * numpy.ones(nx*ny) # Pa

# # Boundary Values

# # xneg
# xneg = numpy.column_stack((xx[:,0].flatten(), yy[:,0].flatten()))
# u_x_xneg = u_x[:,0]
# u_y_xneg = u_y[:,0]
# p_xneg = p[:,0]
# trace_strain_xneg = trace_strain[:,0]

# # xpos
# xpos = numpy.column_stack((xx[:,-1].flatten(), yy[:,-1].flatten()))

# u_x_xpos = u_x[:,-1]
# u_y_xpos = u_y[:,-1]
# p_xpos = p[:,-1]
# trace_strain_xpos = trace_strain[:,-1]

# # yneg
# yneg = numpy.column_stack((xx[0,:].flatten(), yy[0,:].flatten()))
# u_x_yneg = u_x[0,:]
# u_y_yneg = u_y[0,:]
# p_yneg = p[0,:]
# trace_strain_yneg = trace_strain[0,:]

# # ypos
# ypos = numpy.column_stack((xx[-1,:].flatten(), yy[-1,:].flatten()))
# u_x_ypos = u_x[-1,:]
# u_y_ypos = u_y[-1,:]
# p_ypos = p[-1,:]
# trace_strain_ypos = trace_strain[-1,:]



# End of file

