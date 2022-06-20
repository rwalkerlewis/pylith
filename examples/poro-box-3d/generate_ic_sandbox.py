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
t0 = 0.0
dt = 0.001
elapsed = 2.0
tsteps = numpy.int32(elapsed/dt)
t = numpy.linspace(-dt,elapsed,tsteps+2)

def et_disp_x(x, t):
    return x*1E-8*(5*numpy.cos(t*2*numpy.pi) + 2*numpy.cos((t-0.5)*2*numpy.pi) + 1*numpy.cos((t+0.3)*0.5*numpy.pi))

def et_disp_y(y, t):
    return y*1E-8*(7*numpy.cos(t*2*numpy.pi) + 4*numpy.cos((t-0.3)*2*numpy.pi) + 7*numpy.cos((t+0.6)*0.5*numpy.pi))

def et_disp_z(z, t):
    return z*1E-8*(7*numpy.cos((t-0.5)*2*numpy.pi) + 4*numpy.cos((t-0.8)*2*numpy.pi) + 7*numpy.cos((t+0.1)*4*numpy.pi))

coeff_x = (5*numpy.cos(t*2*numpy.pi) + 2*numpy.cos((t-0.5)*2*numpy.pi) + 1*numpy.cos((t+0.3)*0.5*numpy.pi))
coeff_y = (7*numpy.cos(t*2*numpy.pi) + 4*numpy.cos((t-0.3)*2*numpy.pi) + 7*numpy.cos((t+0.6)*0.5*numpy.pi))
coeff_z = (7*numpy.cos((t-0.5)*2*numpy.pi) + 4*numpy.cos((t-0.8)*2*numpy.pi) + 7*numpy.cos((t+0.1)*4*numpy.pi))

tide_coeff = numpy.column_stack((coeff_x, coeff_y, coeff_z))
time_history = tide_coeff / tide_coeff[0, :]

et_x = et_disp_x(x, t.reshape([t.size, 1]))
et_y = et_disp_y(y, t.reshape([t.size, 1]))
et_z = et_disp_z(z, t.reshape([t.size, 1]))

# Boundary Values

# xneg
xneg = numpy.column_stack((xxx[:,0,:].flatten(), yyy[:,0,:].flatten(), zzz[:,0,:].flatten()))
xneg_disp_x = et_disp_x(xneg[:,0], t0)
xneg_disp_y = et_disp_y(xneg[:,1], t0)
xneg_disp_z = et_disp_z(xneg[:,2], t0)

# xpos
xpos = numpy.column_stack((xxx[:,-1,:].flatten(), yyy[:,-1,:].flatten(), zzz[:,-1,:].flatten()))
xpos_disp_x = et_disp_x(xpos[:,0], t0)
xpos_disp_y = et_disp_y(xpos[:,1], t0)
xpos_disp_z = et_disp_z(xpos[:,2], t0)

# yneg
yneg = numpy.column_stack((xxx[0,:,:].flatten(), yyy[0,:,:].flatten(), zzz[0,:,:].flatten()))
yneg_disp_x = et_disp_x(yneg[:,0], t0)
yneg_disp_y = et_disp_y(yneg[:,1], t0)
yneg_disp_z = et_disp_z(yneg[:,2], t0)

# ypos
ypos = numpy.column_stack((xxx[-1,:,:].flatten(), yyy[-1,:,:].flatten(), zzz[1,:,:].flatten()))
ypos_disp_x = et_disp_x(ypos[:,0], t0)
ypos_disp_y = et_disp_y(ypos[:,1], t0)
ypos_disp_z = et_disp_z(ypos[:,2], t0)

# zneg
zneg = numpy.column_stack((xxx[:,:,0].flatten(), yyy[:,:,0].flatten(), zzz[:,:,0].flatten()))
zneg_disp_x = et_disp_x(zneg[:,0], t0)
zneg_disp_y = et_disp_y(zneg[:,1], t0)
zneg_disp_z = et_disp_z(zneg[:,2], t0)

# zpos
zpos = numpy.column_stack((xxx[:,:,-1].flatten(), yyy[:,:,-1].flatten(), zzz[:,:,-1].flatten()))
zpos_disp_x = et_disp_x(zpos[:,0], t0)
zpos_disp_y = et_disp_y(zpos[:,1], t0)
zpos_disp_z = et_disp_z(zpos[:,2], t0)

# Generate figure for strain vector

ee_x = (5*numpy.cos(t*2*numpy.pi) + 2*numpy.cos((t-0.5)*2*numpy.pi) + 1*numpy.cos((t+0.3)*0.5*numpy.pi))
ee_y = (7*numpy.cos(t*2*numpy.pi) + 4*numpy.cos((t-0.3)*2*numpy.pi) + 7*numpy.cos((t+0.6)*0.5*numpy.pi))
ee_z = (7*numpy.cos((t-0.5)*2*numpy.pi) + 4*numpy.cos((t-0.8)*2*numpy.pi) + 7*numpy.cos((t+0.1)*4*numpy.pi))

ee_v = numpy.column_stack((ee_x, ee_y, ee_z))
# End of file

fig, ax = plt.subplot()