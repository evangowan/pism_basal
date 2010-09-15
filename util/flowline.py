#!/usr/bin/env python
# update the history attribute

try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

import numpy as np

from sys import argv

filename = argv[1]

def copy_dim(nc1, nc2, name, direction):
    """Copy a dimension from nc1 to nc2."""
    dim1 = nc1.dimensions[name]
    dim2 = nc2.createDimension(name, len(dim1))

def copy_attributes(var1, var2):
    """Copy attributes of var1 to var2."""
    for each in var1.ncattrs():
        var2.setncattr(each, var1.getncattr(each))


def process(input, output, direction, action="collapse"):
    """Process the file 'input', expanding or collapsing data according to
    'action' and 'direction'. Saves result in 'output'."""
    nc = CDF(input)

    out = CDF(output, 'w', format="NETCDF3_CLASSIC")
    copy_attributes(nc, out)

    for name in nc.dimensions.keys():
        copy_dim(nc, out, name, direction)

    if action == "collapse":
        for name in nc.variables.keys():
            collapse_var(nc, out, name, direction)
    else:
        for name in nc.variables.keys():
            expand_var(nc, out, name, direction)

    out.close()

def collapse_var(nc, out, name, direction):
    """Saves a collapsed (according to 'direction')
    copy of a variable 'name' in 'nc' to 'out'."""
    var1 = nc.variables[name]
    N = (len(nc.dimensions[direction]) - 1) / 2

    print "Processing", name
    dims = var1.dimensions
    if len(dims) > 1:                   # only collapse spatial fields
        dims = filter(lambda(x): x != direction, dims)
        
    var2 = out.createVariable(name, var1.dtype, dimensions=dims)
    copy_attributes(var1, var2)

    # Note: the following depends on the current PISM variable storage order: (t, z, y, x)
    if direction == 'x':
        if var1.ndim == 4:
            var2[:] = var1[:, :, :, N]  # (t, z, y, x)
        elif var1.ndim == 3:
            var2[:] = var1[:, :, N]     # (t, y, x)
        elif var1.ndim == 2:
            var2[:] = var1[:, N]        # (y, x)
        elif var1.ndim == 1:
            var2[:] = var1[:]           # (t), (y) or (x)
    elif direction == 'y':
        if var1.ndim == 4:
            var2[:] = var1[:, :, N, :]
        elif var1.ndim == 3:
            var2[:] = var1[:, N, :]
        elif var1.ndim == 2:
            var2[:] = var1[N, :]
        elif var1.ndim == 1:
            var2[:] = var1[:]

def expand_var(nc, out, name, direction):
    """Saves an expanded (according to 'direction')
    copy of a variable 'name' in 'nc' to 'out'."""
    var1 = nc.variables[name]
    length = len(nc.dimensions[direction])

    print "Processing", name

    # Copy coordinate variables and stop:
    if name in ['t', 'z', 'y', 'x', 'zb']:
        var2 = out.createVariable(name, var1.dtype, (name,))
        var2[:] = var1[:]
        return

    dims = var1.dimensions
    if   len(dims) == 1:
        dims = ('y', 'x')
    elif len(dims) == 2:
        dims = ('t', 'y', 'x')
    elif len(dims) == 3:
        if name == "litho_temp":        # litho_temp is the only variable depending on 'zb'.
            dims = ('t', 'zb', 'y', 'x')
        else:
            dims = ('t', 'z', 'y', 'x')

    var2 = out.createVariable(name, var1.dtype, dims)
    copy_attributes(var1, var2)

    # Note: the following depends on the current PISM variable storage order: (t, z, y, x)
    if direction == 'x':
        if var1.ndim == 3:
            for j in range(length):
                var2[:,:,:,j] = var1[:,:,:]
        elif var1.ndim == 2:
            for j in range(length):
                var2[:,:,j] = var1[:,:]
        elif var1.ndim == 1:
            for j in range(length):
                var2[:,j] = var1[:]
    elif direction == 'y':
        if var1.ndim == 3:
            for j in range(length):
                var2[:,:,j,:] = var1[:,:,:]
        elif var1.ndim == 2:
            for j in range(length):
                var2[:,j,:] = var1[:,:]
        elif var1.ndim == 1:
            for j in range(length):
                var2[j,:] = var1[:]

process(filename, "foo.nc", 'y')
process("foo.nc", "bar.nc", 'y', action="expand")
