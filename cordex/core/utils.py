import tempfile
import warnings

import xarray as xr

from . import cf


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


def pole(ds):
    """Returns rotated pole longitude and latitude"""
    pole_lon = ds[cf.DEFAULT_MAPPING_NCVAR].grid_north_pole_longitude
    pole_lat = ds[cf.DEFAULT_MAPPING_NCVAR].grid_north_pole_latitude
    return pole_lon, pole_lat


def pole_crs(ds):
    """Return a cartropy RotatedPole instance"""
    from cartopy.crs import RotatedPole

    return RotatedPole(*pole(ds))


def _map_crs(x_stack, y_stack, src_crs, trg_crs=None):
    """coordinate transformation of longitude and latitude"""

    from cartopy import crs as ccrs

    if trg_crs is None:
        trg_crs = ccrs.PlateCarree()
    result = trg_crs.transform_points(src_crs, x_stack, y_stack)
    return result[:, :, 0], result[:, :, 1]


# wrapper function for xarray.apply_ufunc
def map_crs(x, y, src_crs, trg_crs=None):
    """coordinate transformation using cartopy

    Transforms the coordinates x, y from the source crs
    into the target crs using cartopy.

    Parameters
    ----------
    x : float array like
        x coordinate of source crs.
    y : float array like
        y coordinate of source crs.
    src_crs : cartopy.crs
        Source coordinate reference system in which x and y
        are defined.
    trg_crs : cartopy.crs
        Target coordinate reference system into which x and y
        should be transformed. If `None`, `PlateCarree` is used.

    Returns
    -------
    x_map : xr.DataArray
        Projected x coordinate.
    y_map : xr.DataArray
        Projected y coordinate.

    """
    warnings.warn("output shape has changed to apply to COARDS conventions")
    y_stack, x_stack = xr.broadcast(y, x)
    input_core_dims = 2 * [list(x_stack.dims)] + [[], []]
    output_core_dims = 2 * [list(x_stack.dims)]
    result = xr.apply_ufunc(
        _map_crs,  # first the function
        x_stack,  # now arguments in the order expected by 'interp1_np'
        y_stack,
        src_crs,
        trg_crs,
        input_core_dims=input_core_dims,  # list with one entry per arg
        # [["rlat", "rlon"], ["rlat", "rlon"]],
        output_core_dims=output_core_dims
        # exclude_dims=set(("lat",)),  # dimensions allowed to change size. Must be set!
    )
    result[0].name = "x_map"
    result[1].name = "y_map"
    return result


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
