# -*- coding: utf-8 -*-
# flake8: noqa
"""Domain module

This module defines preconfigured CORDEX domain from csv tables. The module
also contains some tools to create a domain dataset from a csv tables or simply
from grid information.

Example:

    To get a list of available implementations, create cordex domains, write
    them to netcdf with some dummy data, you can use ,e.g.,::

        from cordex import domain as dm

        eur11 = dm.cordex_domain('EUR-11')

"""

import numpy as np
import pandas as pd
import xarray as xr

from ..tables import domains
from . import cf
from . import utils


def domain_names(table_name=None):
    """Returns a list of short names of all availabe Cordex domains

    Parameters
    ----------
    table_name:
        Only return domain names from this table.

    Returns
    -------
    domain names : list
        List of available cordex domains.

    """
    if table_name:
        return list(domains.tables[table_name].index)
    else:
        return list(domains.table.index)


def cordex_domain(
    short_name,
    dummy=False,
    add_vertices=False,
    tables=None,
    attrs=None,
    mapping_name=None,
    bounds=None,
):
    """Creates an xarray dataset containg the domain grid definitions.

    Parameters
    ----------
    short_name:
        Name of the Cordex Domain.
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.
    add_vertices : bool
        Add grid boundaries in the gloabl coordinates (lon_vertices and lat_vertices).
    tables: dataframe or list of dataframes, default: cordex_tables
        Tables from which to look up the grid information. Index in the table
        should be the short name of the domain, e.g., `EUR-11`. If no table is
        provided, all standard tables will be searched.
    attrs: str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name: str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    bounds: str
        Add spatial bounds.

    Returns
    -------
    Dataset : xarray.core.Dataset
        Dataset containing the coordinates.

    Example
    -------

    To create a cordex rotated pole domain dataset, you can use ,e.g.,::

        import cordex as cx

        eur11 = cx.cordex_domain('EUR-11')

    """
    if attrs is None:
        attrs = {}
    if tables is None:
        tables = domains.table
    if isinstance(tables, list):
        config = pd.concat(tables).loc[short_name]
    else:
        config = tables.loc[short_name]
    return create_dataset(
        **config,
        name=short_name,
        dummy=dummy,
        add_vertices=add_vertices,
        attrs=attrs,
        mapping_name=mapping_name,
        bounds=bounds
    )


def create_dataset(
    nlon,
    nlat,
    dlon,
    dlat,
    ll_lon,
    ll_lat,
    pollon,
    pollat,
    name=None,
    dummy=False,
    add_vertices=False,
    attrs=None,
    mapping_name=None,
    bounds=None,
    **kwargs
):
    """Create domain dataset from grid information.

    Parameters
    ----------
    nlon : int
        longitudal number of grid boxes
    nlat : int
        latitudal number of grid boxes
    dlon : float
        longitudal resolution (degrees)
    dlat : float
        latitudal resolution (degrees)
    ll_lon : float
        lower left rotated longitude (degrees)
    ll_lat : float
        lower left rotated latitude (degrees)
    pollon : float
        pol longitude (degrees)
    pollat : float
        pol latitude (degrees)
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.
    add_vertices : bool
        Add grid boundaries in the gloabl coordinates (lon_vertices and lat_vertices).
    attrs: str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name: str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    """
    if attrs == "CORDEX":
        attrs = cf.DEFAULT_CORDEX_ATTRS
    elif attrs is None:
        attrs = {}
    if name:
        attrs["CORDEX_domain"] = name
    rlon, rlat = _init_grid(nlon, nlat, dlon, dlat, ll_lon, ll_lat)
    lon, lat = rotated_coord_transform(*_stack(rlon, rlat), pollon, pollat)
    pole = _grid_mapping(pollon, pollat)
    return _get_dataset(
        rlon,
        rlat,
        lon,
        lat,
        pole,
        add_vertices=add_vertices,
        dummy=dummy,
        mapping_name=mapping_name,
        attrs=attrs,
    )


def domain_info(short_name, tables=None):
    """Returns a dictionary containg the domain grid definitions.

    Returns a dictionary with grid information according to the
    Cordex archive specifications.

    See https://is-enes-data.github.io/cordex_archive_specifications.pdf

    Parameters
    ----------
    short_name:
        Name of the Cordex Domain.

    Returns
    -------
    domain info : dict
        Dictionary containing the grid information.

    """
    if tables is None:
        tables = domains.table
    if isinstance(tables, list):
        config = pd.concat(tables).loc[short_name]
    else:
        config = tables.loc[short_name]
    return {**{"short_name": short_name}, **dict(**config)}


def _get_dataset(
    rlon,
    rlat,
    lon,
    lat,
    pole,
    add_vertices=False,
    dummy=None,
    mapping_name=None,
    attrs=None,
):
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    data_vars = {mapping_name: pole}

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            rlon=(cf.RLON_NAME, rlon),
            rlat=([cf.RLAT_NAME], rlat),
            lon=([cf.RLAT_NAME, cf.RLON_NAME], lon),
            lat=([cf.RLAT_NAME, cf.RLON_NAME], lat),
        ),
        attrs=attrs,
    )

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cf.coords[key]

    if add_vertices is True:
        from cartopy import crs as ccrs

        pole = (
            ds[mapping_name].grid_north_pole_longitude,
            ds[mapping_name].grid_north_pole_latitude,
        )
        v = vertices(ds.rlon, ds.rlat, ccrs.RotatedPole(*pole))
        ds = xr.merge([ds, v])
        ds[cf.LON_NAME].attrs["bounds"] = cf.LON_BOUNDS
        ds[cf.LAT_NAME].attrs["bounds"] = cf.LAT_BOUNDS

    if dummy:
        if dummy is True:
            dummy_name = "dummy"
        else:
            dummy_name = dummy
        dummy = xr.DataArray(
            data=np.zeros((ds.rlat.size, ds.rlon.size)),
            coords=(ds.rlat, ds.rlon),
        )
        dummy.attrs = {"grid_mapping": mapping_name, "coordinates": "lat lon"}
        ds[dummy_name] = dummy
        if dummy_name == "topo":
            # use cdo to create dummy topography data.
            from cdo import Cdo

            tmp = utils.get_tempfile()
            ds.to_netcdf(tmp)
            topo = Cdo().topo(tmp, returnXDataset=True)["topo"]
            ds[dummy_name] = xr.DataArray(
                data=topo.values,
                coords=(ds.rlat, ds.rlon),
            )
            ds[dummy_name].attrs.update(topo.attrs)
    return ds


def _grid_mapping(pollon, pollat, mapping_name=None):
    """creates a grid mapping DataArray object"""
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    da = xr.DataArray(np.zeros((), dtype=np.int32))
    attrs = cf.mapping.copy()
    attrs["grid_north_pole_longitude"] = pollon
    attrs["grid_north_pole_latitude"] = pollat
    da.attrs = attrs
    return da


def _init_grid(nlon, nlat, dlon, dlat, ll_lon, ll_lat):
    """create coordinate arrays from lower left longitude and latitude"""
    rlon = np.array([ll_lon + i * dlon for i in range(0, nlon)], dtype=np.float64)
    rlat = np.array([ll_lat + i * dlat for i in range(0, nlat)], dtype=np.float64)
    return rlon, rlat


def _stack(x, y):
    """Stack 1d arrays into 2d fields."""
    tmp_x = np.array(x).squeeze()
    tmp_y = np.array(y).squeeze()
    return np.vstack(len(tmp_y) * (tmp_x,)), np.hstack(
        len(tmp_x) * (tmp_y[:, np.newaxis],)
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


def _map_crs(lon, lat, src_crs, trg_crs=None):
    """coordinate transformation of longitude and latitude

    Transforms the coordinates lat, lon from the transform crs
    into the projection crs using cartopy.crs.

    Parameters
    ----------
    lon : float array like
        Longitude coordinate.
    lat : float array like
        Latitude coordinate.
    src_crs : cartopy.crs
        Source coordinate reference system into which lat and lon
        should are defined.
    trg_crs : cartopy.crs
        Target coordinate reference system in which lat and lon
        should be transformed. If `None`, `PlateCarree` is used.

    Returns
    -------
    lon : array like
        Projected longitude coordinate.
    lat : array like
        Projected latitude coordinate.

    """

    from cartopy import crs as ccrs

    if trg_crs is None:
        trg_crs = ccrs.PlateCarree()
    # latlon = ccrs.PlateCarree()
    # lon_stack, lat_stack = xr.broadcast(lon, lat)
    lon_stack = np.broadcast_to(lon, (lat.shape[0], lon.shape[0])).T
    lat_stack = np.broadcast_to(lat, (lon.shape[0], lat.shape[0]))
    result = trg_crs.transform_points(src_crs, lon_stack, lat_stack)
    return result[:, :, 0], result[:, :, 1]


# wrapper function for xarray.apply_ufunc
def map_crs(lon, lat, src_crs, trg_crs=None):
    """coordinate transformation of longitude and latitude

    Transforms the coordinates lat, lon from the transform crs
    into the projection crs using cartopy.crs.

    Parameters
    ----------
    lon : float array like
        Longitude coordinate.
    lat : float array like
        Latitude coordinate.
    src_crs : cartopy.crs
        Source coordinate reference system into which lat and lon
        should are defined.
    trg_crs : cartopy.crs
        Target coordinate reference system in which lat and lon
        should be transformed. If `None`, `PlateCarree` is used.

    Returns
    -------
    lon : xr.DataArray
        Projected longitude coordinate with dims (lat, lon).
    lat : xr.DataArray
        Projected latitude coordinate with dims (lat, lon).

    """
    input_core_dims = [[lon.dims[0]], [lat.dims[0]]] + [[], []]
    output_core_dims = 2 * [[lon.dims[0], lat.dims[0]]]
    result = xr.apply_ufunc(
        _map_crs,  # first the function
        lon,  # now arguments in the order expected by 'interp1_np'
        lat,
        src_crs,
        trg_crs,
        input_core_dims=input_core_dims,  # list with one entry per arg
        output_core_dims=output_core_dims  # [["rlat", "rlon"], ["rlat", "rlon"]],
        # exclude_dims=set(("lat",)),  # dimensions allowed to change size. Must be set!
    )
    result[0].name = cf.LON_NAME
    result[1].name = cf.LAT_NAME
    return result


def _dcoord(coord, include="left"):
    dcoord = coord.values[1:] - coord.values[:-1]
    if include in ["left", "both"]:
        dcoord = np.insert(dcoord, 0, dcoord[0])
    if include in ["right", "both"]:
        dcoord = np.insert(dcoord, dcoord.size, dcoord[-1])
    return dcoord


def _bounds(coord, include="left"):
    dc = _dcoord(coord, include)
    left = coord - 0.5 * dc
    right = coord + 0.5 * dc
    left.name = "left"
    right.name = "right"
    return xr.merge([left, right])


def bounds_coordinates(ds, coords):
    """Adds coordinate bounds as coordinates to the dataset."""
    if isinstance(coords, str):
        coords = (coords,)
    bounds = []
    for coord in coords:
        lr = _bounds(ds.coords[coord])
        name = ds.coords[coord].name + "_b"
        bounds.append(
            xr.DataArray(np.append(lr.left, lr.right[-1]), dims=name, name=name)
        )
    return ds.assign_coords({b.name: b for b in bounds})


def bounds(coords):
    if isinstance(coords, xr.DataArray):
        coords = (coords,)
    return xr.merge([_bounds(coord).left for coord in coords])
    # return xr.merge([_bounds(coord).left.rename({"left": coord.name+"bounds"}) for coord in coords])


def vertices(rlon, rlat, src_crs, trg_crs=None):
    """Compute lon and lat vertices.

    Transformation of rlon vertices and rlat vertices
    into the target crs according to
    https://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries

    Parameters
    ----------
    rlon : xr.DataArray
        Longitude in rotated pole grid.
    rlat : xr.DataArray
        Latitude in rotated pole grid.

    Returns
    -------
    vertices : xr.Dataset
        lon_vertices and lat_vertices in target crs.

    """
    rlon_bounds = _bounds(rlon)
    rlat_bounds = _bounds(rlat)
    # maps each vertex to lat lon coordinates
    # order is counterclockwise starting from lower left vertex
    v1 = map_crs(rlon_bounds.left, rlat_bounds.left, src_crs, trg_crs)
    v2 = map_crs(rlon_bounds.right, rlat_bounds.left, src_crs, trg_crs)
    v3 = map_crs(rlon_bounds.right, rlat_bounds.right, src_crs, trg_crs)
    v4 = map_crs(rlon_bounds.left, rlat_bounds.right, src_crs, trg_crs)
    lon_vertices = xr.concat(
        [v1[0], v2[0], v3[0], v4[0]], dim=cf.BOUNDS_DIM
    ).transpose()
    #    ..., "vertices"
    # )
    lat_vertices = xr.concat(
        [v1[1], v2[1], v3[1], v4[1]], dim=cf.BOUNDS_DIM
    ).transpose()
    #    ..., "vertices"
    # )
    lon_vertices.name = cf.LON_BOUNDS
    lat_vertices.name = cf.LAT_BOUNDS
    lon_vertices.attrs = cf.coords[cf.LON_BOUNDS]
    lat_vertices.attrs = cf.coords[cf.LAT_BOUNDS]
    return xr.merge([lat_vertices, lon_vertices])


def lon_lat_bounds_coordinates(ds, src_crs=None, trg_crs=None):
    """Compute lon and lat bounds coordinates.

    Transformation of rlon vertices and rlat vertices
    into the target crs according to
    https://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries

    Parameters
    ----------
    rlon : xr.DataArray
        Longitude in rotated pole grid.
    rlat : xr.DataArray
        Latitude in rotated pole grid.

    Returns
    -------
    vertices : xr.Dataset
        lon_vertices and lat_vertices in target crs.

    """
    if src_crs is None:
        src_crs = get_pole_crs(ds)
    rot_bounds = bounds_coordinates(ds, ("rlon", "rlat"))
    lon_bounds, lat_bounds = map_crs(
        rot_bounds.coords["rlon_b"], rot_bounds.coords["rlat_b"], src_crs, trg_crs
    )
    return ds.assign_coords(lon_b=lon_bounds.transpose(), lat_b=lat_bounds.transpose())


def get_pole_crs(ds, grid_mapping_name="rotated_latitude_longitude"):
    """Create cartopy crs from rotated pole dataset."""
    from cartopy import crs as ccrs

    pole = (
        ds[grid_mapping_name].grid_north_pole_longitude,
        ds[grid_mapping_name].grid_north_pole_latitude,
    )
    return ccrs.RotatedPole(*pole)


def extend_coordinate(coord, n=1):
    data = coord
    # insert = np.array([data[0] - n*(data[1]-data[0]) + i*(data[1]-data[0]) for i in range(0, n)])
    insert = np.array([data[0] - (n - i) * (data[1] - data[0]) for i in range(0, n)])
    return insert
    data_extended = np.insert(data, 0, data[0] - (data[1] - data[0]))
    # return data_extended
    # data_extended = np.append(data_extended, data[-1]+(data[-1]-data[-2]))
    return xr.DataArray(data_extended, dims=coord.dims, name=coord.name)
