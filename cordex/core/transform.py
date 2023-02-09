from warnings import warn

import xarray as xr
from pyproj import CRS, Transformer


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
    warn(
        "map_crs is deprecated, please use transform instead",
        DeprecationWarning,
        stacklevel=2,
    )
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


def _transform(x, y, src_crs, trg_crs):
    """helper function for transforming coordinates"""
    transformer = Transformer.from_crs(src_crs, trg_crs)
    yt, xt = transformer.transform(x, y)
    return xt, yt


def transform(x, y, src_crs, trg_crs=None):
    """Coordinate transformation using pyproj.

    Transforms the coordinates x, y from the source crs
    into a target crs using pyproj.

    Parameters
    ----------
    x : DataArray
        Longitude coordinate.
    y : DataArray
        Latitude coordinate.
    src_crs : pyproj.crs
        Source coordinate reference system into which x and y are defined.
    trg_crs : pyproj.crs
        Target coordinate reference system into which x and y
        should be transformed. If `None`, `EPSG:4326` is used.

    Returns
    -------
    xt : DataArray
        Transformed x coordinate.
    yt : DataArray
        Transformed y coordinate.

    """
    if trg_crs is None:
        # default target crs
        trg_crs = CRS("EPSG:4326")
    y_stack, x_stack = xr.broadcast(y, x)
    input_core_dims = [x_stack.dims, y_stack.dims] + [[], []]
    output_core_dims = [x_stack.dims, y_stack.dims]
    xt, yt = xr.apply_ufunc(
        _transform,
        x_stack,
        y_stack,
        src_crs,
        trg_crs,
        input_core_dims=input_core_dims,
        output_core_dims=output_core_dims,
    )
    xt.name = "xt"
    yt.name = "yt"
    return xt, yt
