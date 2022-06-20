#!/usr/bin/env nemesis
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

from time import time
import numpy 

# Dimensions
DOMAIN_X = 20.0
DOMAIN_Y = 20.0
DOMAIN_Z = 20.0

# Discretization
dx = 0.25
dy = 0.25 
dz = 0.25

nx = numpy.int64(DOMAIN_X / dx)
ny = numpy.int64(DOMAIN_Y / dy)
nz = numpy.int64(DOMAIN_Z / dz)

# Grid
x = numpy.linspace(-DOMAIN_X/2,DOMAIN_X/2,nx+1)
y = numpy.linspace(-DOMAIN_Y/2,DOMAIN_Y/2,ny+1)
z = numpy.linspace(-DOMAIN_Z/2,DOMAIN_Z/2,nz+1)

xxx, yyy, zzz = numpy.meshgrid(x,y,z)
xyz = numpy.column_stack((xxx.flatten(), yyy.flatten(), zzz.flatten()))

# Time in decimal days
dt = 0.01
elapsed = 3.0
tsteps = numpy.int32(elapsed/dt)
t = numpy.linspace(-dt,elapsed,tsteps+2)

# Earth Tide Values

# Time in decimal days
t0 = 0.0
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

class GenerateDB(object):
    def run(self):
        """Generate the databases.
        """
        # displacement
        self.xneg_disp_x()
        self.xneg_disp_y()
        self.xneg_disp_z()                

        self.xpos_disp_x()
        self.xpos_disp_y()
        self.xpos_disp_z()                

        self.yneg_disp_x()
        self.yneg_disp_y()
        self.yneg_disp_z()

        self.ypos_disp_x()
        self.ypos_disp_y()
        self.ypos_disp_z()

        self.zneg_disp_x()
        self.zneg_disp_y()
        self.zneg_disp_z()                        

        self.zpos_disp_x()
        self.zpos_disp_y()
        self.zpos_disp_z()

        # time history
        self.time_history_x()
        self.time_history_y()
        self.time_history_z()        

        return

# ========================================================================

    def xneg_disp_x(self):
        """Generate the database for xneg displacement, x component.
        """
        # X Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(xneg_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xneg_disp_z.size)}

        data = {"points": xneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xneg_disp_x.spatialdb")
        io.write(data)        
        return

    def xneg_disp_y(self):
        """Generate the database for xneg displacement, y component.
        """
        # X Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(xneg_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xneg_disp_z.size)}

        data = {"points": xneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xneg_disp_y.spatialdb")
        io.write(data)        
        return

    def xneg_disp_z(self):
        """Generate the database for xneg displacement, z component.
        """
        # X Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                               "units": "m",
                               "data": numpy.zeros(xneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(xneg_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xneg_disp_z.size)}                        

        data = {"points": xneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xneg_disp_z.spatialdb")
        io.write(data)        
        return

    def xpos_disp_x(self):
        """Generate the database for xpos displacement, x component.
        """
        # X Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(xpos_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xpos_disp_z.size)}

        data = {"points": xpos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xpos_disp_x.spatialdb")
        io.write(data)        
        return

    def xpos_disp_y(self):
        """Generate the database for xpos displacement, y component.
        """
        # X Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(xpos_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xpos_disp_z.size)}

        data = {"points": xpos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xpos_disp_y.spatialdb")
        io.write(data)        
        return

    def xpos_disp_z(self):
        """Generate the database for xpos displacement, z component.
        """
        # X Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                               "units": "m",
                               "data": numpy.zeros(xpos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(xpos_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(xpos_disp_z.size)}

        data = {"points": xpos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_xpos_disp_z.spatialdb")
        io.write(data)        
        return

# ========================================================================

    def yneg_disp_x(self):
        """Generate the database for yneg displacement, x component.
        """
        # Y Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(yneg_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                               "units": "m",
                               "data": numpy.zeros(yneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                               "units": "m",
                               "data": numpy.zeros(yneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(yneg_disp_z.size)}

        data = {"points": yneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_yneg_disp_x.spatialdb")
        io.write(data)        
        return

    def yneg_disp_y(self):
        """Generate the database for yneg displacement, y component.
        """
        # Y Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(yneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(yneg_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(yneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(yneg_disp_z.size)}

        data = {"points": yneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_yneg_disp_y.spatialdb")
        io.write(data)        
        return

    def yneg_disp_z(self):
        """Generate the database for yneg displacement, z component.
        """
        # X Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(yneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(yneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(yneg_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(yneg_disp_z.size)}

        data = {"points": yneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start
                ]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_yneg_disp_z.spatialdb")
        io.write(data)        
        return

    def ypos_disp_x(self):
        """Generate the database for ypos displacement, x component.
        """
        # Y Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(ypos_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(ypos_disp_z.size)}

        data = {"points": ypos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_ypos_disp_x.spatialdb")
        io.write(data)        
        return

    def ypos_disp_y(self):
        """Generate the database for ypos displacement, y component.
        """
        # Y Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(ypos_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(ypos_disp_z.size)}

        data = {"points": ypos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_ypos_disp_y.spatialdb")
        io.write(data)        
        return

    def ypos_disp_z(self):
        """Generate the database for ypos displacement, z component.
        """
        # Y Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(ypos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(ypos_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(ypos_disp_z.size)}

        data = {"points": ypos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_ypos_disp_z.spatialdb")
        io.write(data)        
        return

# ========================================================================

    def zneg_disp_x(self):
        """Generate the database for zneg displacement, x component.
        """
        # Z Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(zneg_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zneg_disp_z.size)}

        data = {"points": zneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zneg_disp_x.spatialdb")
        io.write(data)        
        return

    def zneg_disp_y(self):
        """Generate the database for zneg displacement, y component.
        """
        # Z Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(zneg_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zneg_disp_z.size)}

        data = {"points": zneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zneg_disp_y.spatialdb")
        io.write(data)        
        return

    def zneg_disp_z(self):
        """Generate the database for zneg displacement, z component.
        """
        # Z Neg Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(zneg_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(zneg_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zneg_disp_z.size)}

        data = {"points": zneg,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zneg_disp_z.spatialdb")
        io.write(data)        
        return

    def zpos_disp_x(self):
        """Generate the database for zpos displacement, x component.
        """
        # Z Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.ravel(zpos_disp_x)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zpos_disp_z.size)}

        data = {"points": zpos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zpos_disp_x.spatialdb")
        io.write(data)        
        return

    def zpos_disp_y(self):
        """Generate the database for zpos displacement, y component.
        """
        # Z Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.ravel(zpos_disp_y)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_z.size)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zpos_disp_z.size)}

        data = {"points": zpos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zpos_disp_y.spatialdb")
        io.write(data)        
        return

    def zpos_disp_z(self):
        """Generate the database for zpos displacement, z component.
        """
        # Z Pos Displacement
        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.units = 'meter'
        cs.inventory.spaceDim = 3
        cs._configure()
        
        displacement_x = {"name": "time_history_amplitude_x",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_x.size)}

        displacement_y = {"name": "time_history_amplitude_y",
                        "units": "m",
                        "data": numpy.zeros(zpos_disp_y.size)}

        displacement_z = {"name": "time_history_amplitude_z",
                        "units": "m",
                        "data": numpy.ravel(zpos_disp_z)}

        time_history_start = {"name": "time_history_start_time",
                               "units": "second",
                               "data": numpy.zeros(zpos_disp_z.size)}

        data = {"points": ypos,                
                "coordsys": cs,
                "data_dim": 2,
                "values": [displacement_x, displacement_y, displacement_z, time_history_start]}
        from spatialdata.spatialdb.SimpleIOAscii import createWriter
        io = createWriter("bc_zpos_disp_z.spatialdb")
        io.write(data)        
        return

    def time_history_x(self):
        """Generate the time history database.
        """
        # tmax = 1
        # nsteps = 2001
        # time = numpy.zeros(nsteps+1)
        # time[:-1] = numpy.linspace(0,tmax,nsteps)
        # time[-1] = 999.000
        # amplitude = time.copy()
        from spatialdata.spatialdb.TimeHistoryIO import write
        write(time=t, amplitude=time_history[:,0], units="day", filename="time_history_x.timedb")        
        return

    def time_history_y(self):
        """Generate the time history database.
        """
        # tmax = 1
        # nsteps = 2001
        # time = numpy.zeros(nsteps+1)
        # time[:-1] = numpy.linspace(0,tmax,nsteps)
        # time[-1] = 999.000
        # amplitude = time.copy()
        from spatialdata.spatialdb.TimeHistoryIO import write
        write(time=t, amplitude=time_history[:,1], units="day", filename="time_history_y.timedb")        
        return

    def time_history_z(self):
        """Generate the time history database.
        """
        # tmax = 1
        # nsteps = 2001
        # time = numpy.zeros(nsteps+1)
        # time[:-1] = numpy.linspace(0,tmax,nsteps)
        # time[-1] = 999.000
        # amplitude = time.copy()
        from spatialdata.spatialdb.TimeHistoryIO import write
        write(time=t, amplitude=time_history[:,2], units="day", filename="time_history_z.timedb")        
        return

# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file

