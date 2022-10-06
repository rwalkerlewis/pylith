#!/usr/bin/env nemesis

import h5py
import matplotlib.pyplot as pyplot

#fullpath = "/Users/baagaard/scratch/build/clang-13.0/cig/pylith-debug/tests/manual/petsc/sliderblock/sliderblock.h5"
fullpath = "/opt/pylith/build/debug/pylith/tests/manual/petsc/sliderblock.h5"

h5 = h5py.File(fullpath, "r")
t = h5["time"][:].squeeze()
soln = h5["solution"][:]
u = soln[:, 0:4, :]
v = soln[:, 4:8, :]
l = soln[:, -1, :]

x = [0, 1, 2, 3]


fig, axes = pyplot.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)
for i in range(4):
    axes.plot(t,x[i] + u[:, i, :].squeeze(), label=r'$\mathregular{u_{s1}}$'.replace('s1', str(i)))
axes.set_ylabel('Distance Units')
axes.set_xlabel('Time, sec')
axes.legend()
fig.tight_layout()
fig.savefig('displacement.png',dpi = 300)
pyplot.show()


fig, axes = pyplot.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)
axes.plot(t,(u[:, 2, :]-u[:, 1, :]).squeeze(), label=r'$\mathregular{u_{2} - u_{1}}$')
axes.set_ylabel('Distance Units')
axes.set_xlabel('Time, sec')
axes.legend()
fig.tight_layout()
fig.savefig('u_2_minus_u_1.png',dpi = 300)
pyplot.show()

fig, axes = pyplot.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)
axes.plot(t,(v[:, 2, :]-v[:, 1, :]).squeeze(), label=r'$\mathregular{v_{2} - v_{1}}$')
axes.set_ylabel('Velocity Units')
axes.set_xlabel('Time, sec')
axes.legend()
fig.tight_layout()
fig.savefig('v_2_minus_v_1.png',dpi = 300)
pyplot.show()

fig, axes = pyplot.subplots(ncols=1, nrows=1)
fig.set_size_inches(12,10)
axes.plot(t,l.squeeze())
axes.set_ylabel(r'$\lambda$')
axes.set_xlabel('Time, sec')
#axes.legend()
fig.tight_layout()
fig.savefig('lambda.png',dpi = 300)
pyplot.show()
#pyplot.plot(t, (u[:, 2, :]-u[:, 1, :]).squeeze())
#pyplot.show()

#pyplot.plot(t, (v[:, 2, :]-v[:, 1, :]).squeeze())
#pyplot.show()

#pyplot.plot(t, l.squeeze())
#pyplot.show()