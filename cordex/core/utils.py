#! /usr/bin/python
# coding: utf-8

import tempfile


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


#
#
# def copy_dataset(src, varname=None, timestep=None, destination=None):
#    """copy an existing NetCDF dataset on disk"""
#    from netCDF4 import Dataset
#
#    if varname is None:
#        variables = src.variables
#    else:
#        variables = {varname: src.variables[varname]}
#    if destination is None:
#        destination = "copy.nc"
#    dst = Dataset(destination, "w")
#    # copy attributes
#    for name in src.ncattrs():
#        dst.setncattr(name, getattr(src, name))
#    # copy dimensions
#    for name, dimension in src.dimensions.items():
#        if timestep and name == "time":
#            length = 1
#        else:
#            length = len(dimension) if not dimension.isunlimited() else None
#        # dst.createDimension( name, (len(dimension) if not dimension.isunlimited() else None))
#        dst.createDimension(name, length)
#        # copy all file data except for the excluded
#    for name, variable in variables.items():
#        if hasattr(variable, "_FillValue"):
#            fill_value = getattr(variable, "_FillValue")
#        else:
#            fill_value = None
#        var = dst.createVariable(
#            name, variable.dtype, variable.dimensions, fill_value=fill_value
#        )
#        for attr in variable.ncattrs():
#            if attr != "_FillValue":
#                dst.variables[name].setncattr(attr, getattr(variable, attr))
#        if variable.shape:
#            if timestep is None or "time" not in variable.dimensions:
#                dst.variables[name][:] = src.variables[name][:]
#            else:
#                dst.variables[name][0] = src.variables[name][timestep]
#    return dst
