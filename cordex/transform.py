from warnings import warn

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

from . import cf

xr.set_options(keep_attrs=True)


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
        output_core_dims=output_core_dims,
        # exclude_dims=set(("lat",)),  # dimensions allowed to change size. Must be set!
    )

    result[0].name = "x_map"
    result[1].name = "y_map"

    return result


def derotate_vector(u, v, lon=None, lat=None, pollon=None, pollat=None):
    """Derotate vector components from rotated coordinates.

    The function performs a backward transformation of, e.g., velocity components u and v
    from a rotated spherical coordinate system to a geographical system. If only the
    components u and v are provided, it is assumed, they are DataArrays containing
    a rotated latitude longitude grid mapping and lon lat coordinates that are used
    for the transformation.

    Parameters
    ----------
    u : float or DataArray
        u component of vector in rotated coordinate system.
    v : float or DataArray
        v component of vector in rotated coordinate system.
    lon : float or DataArray
        Longitude coordinates in which to transform the vector components.
    lat : float or DataArray
        Latitude coordinates in which to transform the vector components.
    pollon : float
        Longitude of north pole in geographical coordinate system.
    pollat : float
        Latitude of north pole in geographical coordinate system.

    Returns
    -------
    ut : DataArray
        Transformed u vector component.
    vt : DataArray
        Transformed v vector component.

    """

    if lon is None:
        lon = u.cf["longitude"]
    if lat is None:
        lat = v.cf["latitude"]

    if pollon is None:
        pollon = u.cf["grid_mapping"].grid_north_pole_longitude
    if pollat is None:
        pollat = v.cf["grid_mapping"].grid_north_pole_latitude

    rla = np.deg2rad(lon)
    phi = np.deg2rad(lat)
    zrla = rla
    zphi = phi

    zpollat = np.deg2rad(pollat)
    zpollon = np.deg2rad(pollon)
    pollond = pollon
    zsinpol = np.sin(zpollat)
    zcospol = np.cos(zpollat)
    if pollon < 0.0:
        pollond = 360.0 + pollon

    # compute vector abs in rotated coordinates

    zzrla = np.rad2deg(rla)
    zzrla = xr.where(zzrla > 180.0, zzrla - 360.0, zzrla)
    zzrla = np.deg2rad(zzrla)
    zd1 = zzrla - zpollon
    zd2 = zrla - zpollon

    # winkel zbeta berechen (schnittwinkel der breitenkreise)

    zarg1 = -np.sin(zd1) * np.cos(zphi)
    zarg2 = -zsinpol * np.cos(zphi) * np.cos(zd1) + zcospol * np.sin(zphi)

    zarg2 = xr.where(abs(zarg2) < 1.0e-20, 1.0e-20, zarg2)
    zrlas = np.arctan2(zarg1, zarg2)
    zarg = -np.sin(zpollat) * np.sin(zd2) * np.sin(zrlas) - np.cos(zd2) * np.cos(zrlas)
    zarg = zarg.clip(min=-1.0, max=1.0)
    zbeta = abs(np.arccos(zarg))

    zbeta = xr.where(
        (-((np.rad2deg(rla)) - (pollond - 180.0)) < 0.0)
        & (-((np.rad2deg(rla)) - (pollond - 180.0)) >= -180.0),
        -zbeta,
        zbeta,
    )

    x1 = u * np.cos(zbeta) - v * np.sin(zbeta)
    y1 = u * np.sin(zbeta) + v * np.cos(zbeta)

    return x1, y1


def _transform(x, y, src_crs, trg_crs):
    """helper function for transforming coordinates"""
    # always_xy=True
    # https://proj.org/faq.html#why-is-the-axis-ordering-in-proj-not-consistent
    transformer = Transformer.from_crs(src_crs, trg_crs, always_xy=True)
    xt, yt = transformer.transform(x, y)
    return xt, yt


def transform(x, y, src_crs, trg_crs=None):
    """Coordinate transformation using pyproj.

    Transforms the coordinates x, y from the source crs
    into a target crs using pyproj.

    Parameters
    ----------
    x : DataArray
        X coordinate.
    y : DataArray
        Y coordinate.
    src_crs : pyproj.CRS
        Source coordinate reference system in which x and y are defined.
    trg_crs : pyproj.CRS
        Target coordinate reference system into which x and y
        should be transformed. If not supplied, ``EPSG:4326`` is the default.

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
    xt.attrs = {"epsg": trg_crs.to_epsg()}
    yt.attrs = {"epsg": trg_crs.to_epsg()}

    return xt, yt


def transform_coords(ds, src_crs=None, trg_crs=None, trg_dims=None):
    """Transform X and Y coordinates of a Dataset.

    The transformed coordinates will be added to the Dataset.
    This function is usefull to add, e.g., global lon/lat coordinates
    to a rotated pole grid.

    Parameters
    ----------
    ds : Dataset or DataArray
        Dataset with input grid.
    src_crs : pyproj.CRS
        Source coordinate reference system in which X and Y are defined.
        If not supplied, a `grid_mapping` variable should be available
        to define the source CRS.
    trg_crs : pyproj.CRS
        Target coordinate reference system into which x and y
        should be transformed. If not supplied, ``EPSG:4326`` is the default.
    trg_dims: list or set
        Names of the output coordinates.

    Returns
    -------
    ds : Dataset or DataArray
        Dataset with transformed coordinates.

    """
    ds = ds.copy(deep=False)

    if trg_crs is None:
        # default target crs
        trg_crs = CRS("EPSG:4326")
    if trg_dims is None:
        trg_dims = ("xt", "yt")
    if src_crs is None:
        src_crs = CRS.from_cf(ds.cf["grid_mapping"].attrs)
    x, y = ds.cf["X"], ds.cf["Y"]
    xt, yt = transform(x, y, src_crs, trg_crs)

    return ds.assign_coords({trg_dims[0]: xt, trg_dims[1]: yt})


def transform_bounds(
    ds, src_crs=None, trg_crs=None, trg_dims=None, bnds_dim=None, keep_xy_bounds=False
):
    """Transform linear X and Y bounds of a Dataset.

    Transformation of of the bounds of linear X and Y coordinates
    into the target crs according to
    https://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    If the linear X and Y coordinate bounds are not available, they will be
    created and transformed.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing linear X and Y coordinates, e.g., `rlon` and `rlat`
    src_crs : pyproj.CRS
        Source coordinate reference system in which X and Y are defined.
        If not supplied, a `grid_mapping` variable should be available
        to define the source CRS.
    trg_crs : pyproj.CRS
        Target coordinate reference system into which bounds
        should be transformed. If not supplied, ``EPSG:4326`` is the default.
    trg_dims: list or set
        Names of the bounds coordinates.
    bnds_dim: str
        Names of the bounds dimension.

    Returns
    -------
    bounds : xr.Dataset
        Dataset with X and Y bounds in target crs. Probably some 2D coordinates.

    References
    ----------
    Please refer to the CF conventions document : https://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries

    """
    ds = ds.copy(deep=False)

    if src_crs is None:
        src_crs = CRS.from_cf(ds.cf["grid_mapping"].attrs)
    if trg_crs is None:
        # default target crs
        trg_crs = CRS("EPSG:4326")
    if trg_dims is None:
        trg_dims = (
            ds.cf["longitude"].name + "_vertices",
            ds.cf["latitude"].name + "_vertices",
        )
    if bnds_dim is None:
        bnds_dim = cf.BOUNDS_DIM

    bnds = ds.cf.add_bounds(("X", "Y"))
    x_bnds = bnds.cf.get_bounds("X").drop(bnds.cf.bounds["X"])
    y_bnds = bnds.cf.get_bounds("Y").drop(bnds.cf.bounds["Y"])

    # order is counterclockwise starting from lower left vertex
    v1 = transform(x_bnds.isel(bounds=0), y_bnds.isel(bounds=0), src_crs, trg_crs)
    v2 = transform(x_bnds.isel(bounds=1), y_bnds.isel(bounds=0), src_crs, trg_crs)
    v3 = transform(x_bnds.isel(bounds=1), y_bnds.isel(bounds=1), src_crs, trg_crs)
    v4 = transform(x_bnds.isel(bounds=0), y_bnds.isel(bounds=1), src_crs, trg_crs)
    xt_vertices = xr.concat([v1[0], v2[0], v3[0], v4[0]], dim=bnds_dim)  # .transpose()
    #    ..., "vertices"
    # )
    yt_vertices = xr.concat([v1[1], v2[1], v3[1], v4[1]], dim=bnds_dim)  # .transpose()

    xt_vertices.name = "xt_vertices"  # cf.LON_BOUNDS
    yt_vertices.name = "yt_vertices"  # cf.LAT_BOUNDS
    xt_vertices.attrs = cf.coords[cf.LON_BOUNDS]
    yt_vertices.attrs = cf.coords[cf.LAT_BOUNDS]

    bounds = xr.merge([xt_vertices, yt_vertices]).transpose(
        ds.cf["Y"].dims[0], ds.cf["X"].dims[0], bnds_dim
    )

    ds[ds.cf["longitude"].name].attrs["bounds"] = trg_dims[0]
    ds[ds.cf["latitude"].name].attrs["bounds"] = trg_dims[1]

    return ds.assign_coords(
        {
            trg_dims[0]: bounds.xt_vertices.drop_vars(
                (ds.cf["X"].name, ds.cf["Y"].name)
            ),
            trg_dims[1]: bounds.yt_vertices.drop_vars(
                (ds.cf["X"].name, ds.cf["Y"].name)
            ),
        }
    )


def rotated_coord_transform(lon, lat, np_lon, np_lat, direction="rot2geo"):
    """Transforms a coordinate into a rotated grid coordinate and vice versa.

    The coordinates have to be given in degree and will be returned in degree.

    Parameters
    ----------
    lon : float array like
        Longitude coordinate.
    lat : float array like
        Latitude coordinate.
    np_lon : float array like
        Longitude coordinate of the rotated north pole.
    np_lat : float array like
        Latitude coordinate of the rotated north pole.
    direction : str
        Direction of the rotation.
        Options are: 'rot2geo' (default) for a transformation to regular
        coordinates from rotated. 'geo2rot' transforms regular coordinates
        to rotated.

    Returns
    -------
    lon_new : array like
        New longitude coordinate.
    lat_new : array like
        New latitude coordinate.
    """
    warn(
        "rotated_coord_transform is deprecated, please use transform_xy instead",
        DeprecationWarning,
        stacklevel=2,
    )
    # Convert degrees to radians
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    theta = 90.0 - np_lat  # Rotation around y-axis
    phi = np_lon + 180.0  # Rotation around z-axis

    # Convert degrees to radians
    phi = np.deg2rad(phi)
    theta = np.deg2rad(theta)

    # Convert from spherical to cartesian coordinates
    x = np.cos(lon) * np.cos(lat)
    y = np.sin(lon) * np.cos(lat)
    z = np.sin(lat)

    # Regular -> Rotated
    if direction == "geo2rot":
        x_new = (
            np.cos(theta) * np.cos(phi) * x
            + np.cos(theta) * np.sin(phi) * y
            + np.sin(theta) * z
        )
        y_new = -np.sin(phi) * x + np.cos(phi) * y
        z_new = (
            -np.sin(theta) * np.cos(phi) * x
            - np.sin(theta) * np.sin(phi) * y
            + np.cos(theta) * z
        )

    # Rotated -> Regular
    elif direction == "rot2geo":
        phi = -phi
        theta = -theta

        x_new = (
            np.cos(theta) * np.cos(phi) * x
            + np.sin(phi) * y
            + np.sin(theta) * np.cos(phi) * z
        )
        y_new = (
            -np.cos(theta) * np.sin(phi) * x
            + np.cos(phi) * y
            - np.sin(theta) * np.sin(phi) * z
        )
        z_new = -np.sin(theta) * x + np.cos(theta) * z

    # Convert cartesian back to spherical coordinates
    lon_new = np.arctan2(y_new, x_new)
    lat_new = np.arcsin(z_new)

    # Convert radians back to degrees
    lon_new = np.rad2deg(lon_new)
    lat_new = np.rad2deg(lat_new)

    return lon_new, lat_new


def grid_mapping(pollon, pollat, mapping_name=None):
    """creates a grid mapping DataArray object"""
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    da = xr.DataArray(np.zeros((), dtype=cf.grid_mapping_dtype))
    attrs = cf.mapping.copy()
    attrs["grid_north_pole_longitude"] = pollon
    attrs["grid_north_pole_latitude"] = pollat
    da.attrs = attrs
    da.name = mapping_name
    return da
