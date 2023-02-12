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

from warnings import warn

import cf_xarray as cfxr
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS

from . import cf
from .config import nround
from .tables import domains
from .transform import grid_mapping, transform, transform_coords
from .utils import get_tempfile


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
        Domain identifier.
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
    pollon=None,
    pollat=None,
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
        lower left longitude (degrees)
    ll_lat : float
        lower left latitude (degrees)
    pollon : float
        pol longitude (degrees)
    pollat : float
        pol latitude (degrees)
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.
    add_vertices : bool
        Add grid boundaries in the global coordinates (lon_vertices and lat_vertices).
    attrs: str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name: str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    """
    rotated = True
    if attrs == "CORDEX":
        attrs = cf.DEFAULT_CORDEX_ATTRS
    elif attrs is None:
        attrs = {}
    if name:
        attrs["CORDEX_domain"] = name
    if pollon is None or pollat is None:
        rotated = False
    try:
        if np.isnan(pollon) or np.isnan(pollat):
            rotated = False
    except:
        pass

    x, y = _lin_coord(nlon, dlon, ll_lon), _lin_coord(nlat, dlat, ll_lat)

    if rotated is True:
        return _get_rotated_dataset(
            x,
            y,
            #  lon,
            #  lat,
            pollon,
            pollat,
            add_vertices=add_vertices,
            dummy=dummy,
            mapping_name=mapping_name,
            attrs=attrs,
        )
    else:
        return _get_regular_dataset(
            x,
            y,
            add_vertices=add_vertices,
            dummy=dummy,
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


def _get_regular_dataset(
    lon,
    lat,
    add_vertices=False,
    dummy=None,
    attrs=None,
):
    ds = xr.Dataset(
        data_vars=None,
        coords=dict(
            lon=(cf.LON_NAME, lon),
            lat=(cf.LAT_NAME, lat),
        ),
        attrs=attrs,
    )

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cf.coords[key]

    ds.lon.attrs["axis"] = "X"
    ds.lat.attrs["axis"] = "Y"

    if add_vertices is True:
        from cartopy import crs as ccrs

        pole = (
            ds[mapping_name].grid_north_pole_longitude,
            ds[mapping_name].grid_north_pole_latitude,
        )
        # v = vertices(ds, ccrs.RotatedPole(*pole))
        v = vertices(ds.rlon, ds.rlat, ccrs.RotatedPole(*pole))
        ds = xr.merge([ds, v])
        ds[cf.LON_NAME].attrs["bounds"] = cf.LON_BOUNDS
        ds[cf.LAT_NAME].attrs["bounds"] = cf.LAT_BOUNDS

    if dummy is not None:
        ds = _add_dummy(ds, dummy)
    return ds


def _get_rotated_dataset(
    rlon,
    rlat,
    pollon,
    pollat,
    add_vertices=False,
    dummy=None,
    mapping_name=None,
    attrs=None,
):
    mapping = grid_mapping(pollon, pollat, mapping_name)
    # lon, lat = transform(x, y, src_crs=CRS.from_cf(mapping.attrs))

    ds = xr.Dataset(
        data_vars={mapping.name: mapping},
        coords=dict(
            rlon=(cf.RLON_NAME, rlon),
            rlat=([cf.RLAT_NAME], rlat),
        ),
        attrs=attrs,
    )

    lon, lat = transform(ds.rlon, ds.rlat, src_crs=CRS.from_cf(mapping.attrs))
    ds = ds.assign_coords(lon=lon, lat=lat)

    if add_vertices is True:
        v = vertices(ds.rlon, ds.rlat, CRS.from_cf(mapping.attrs))
        ds = xr.merge([ds, v])
        ds[cf.LON_NAME].attrs["bounds"] = cf.LON_BOUNDS
        ds[cf.LAT_NAME].attrs["bounds"] = cf.LAT_BOUNDS

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cf.coords[key]

    if dummy is not None:
        ds = _add_dummy(ds, dummy)
    return ds


def _add_dummy(ds, name=True):
    if name is True:
        name = "dummy"

    attrs = {"coordinates": "lat lon"}
    try:
        attrs["grid_mapping"] = ds.cf["grid_mapping"].name
    except KeyError:
        pass
    # this is required for CDO to understand coordinates
    ds[name] = xr.DataArray(
        data=np.zeros((ds.cf.dims["Y"], ds.cf.dims["X"])),
        coords=(ds.cf["Y"], ds.cf["X"]),
        attrs=attrs,
    )
    if name == "topo":
        # use cdo to create dummy topography data.
        from cdo import Cdo

        tmp = get_tempfile()
        ds.to_netcdf(tmp)
        topo = Cdo().topo(tmp, returnXDataset=True)["topo"]
        ds[name] = xr.DataArray(
            data=topo.values, coords=(ds.cf["Y"], ds.cf["X"]), attrs=attrs
        )
        ds[name].attrs.update(topo.attrs)
    return ds


def _lin_coord(nx, dx, x0, dtype=np.float64):
    """create coordinate arrays from lower left longitude and latitude"""
    return np.round(np.arange(0, nx, dtype=dtype) * dx + x0, nround)


def _stack(x, y):
    """Stack 1d arrays into 2d fields."""
    tmp_x = np.array(x).squeeze()
    tmp_y = np.array(y).squeeze()
    return np.vstack(len(tmp_y) * (tmp_x,)), np.hstack(
        len(tmp_x) * (tmp_y[:, np.newaxis],)
    )


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


def vertices_new(ds, src_crs, trg_crs=None):
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
    warn(
        "Order of vertices has changed since v0.3.2 to CF Conventions, see https://github.com/euro-cordex/py-cordex/issues/34"
    )
    # ds = xr.merge([rlon_bounds, rlat_bounds])
    rlon_bounds = ds.cf.add_bounds("rlon").rlon_bounds.drop(
        "rlon_bounds"
    )  # _bounds(ds.rlon)
    rlat_bounds = ds.cf.add_bounds("rlat").rlat_bounds.drop(
        "rlat_bounds"
    )  # _bounds(ds.rlat)
    # maps each vertex to lat lon coordinates
    # order is counterclockwise starting from lower left vertex
    v1 = transform(
        rlon_bounds.isel(bounds=0), rlat_bounds.isel(bounds=0), src_crs, trg_crs
    )
    v2 = transform(
        rlon_bounds.isel(bounds=1), rlat_bounds.isel(bounds=0), src_crs, trg_crs
    )
    v3 = transform(
        rlon_bounds.isel(bounds=1), rlat_bounds.isel(bounds=1), src_crs, trg_crs
    )
    v4 = transform(
        rlon_bounds.isel(bounds=0), rlat_bounds.isel(bounds=1), src_crs, trg_crs
    )
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
    v1 = transform(rlon_bounds.left, rlat_bounds.left, src_crs, trg_crs)
    v2 = transform(rlon_bounds.right, rlat_bounds.left, src_crs, trg_crs)
    v3 = transform(rlon_bounds.right, rlat_bounds.right, src_crs, trg_crs)
    v4 = transform(rlon_bounds.left, rlat_bounds.right, src_crs, trg_crs)
    lon_vertices = xr.concat(
        [v1[0].T, v2[0].T, v3[0].T, v4[0].T], dim=cf.BOUNDS_DIM
    ).transpose()
    #    ..., "vertices"
    # )
    lat_vertices = xr.concat(
        [v1[1].T, v2[1].T, v3[1].T, v4[1].T], dim=cf.BOUNDS_DIM
    ).transpose()
    #    ..., "vertices"
    # )
    lon_vertices.name = cf.LON_BOUNDS
    lat_vertices.name = cf.LAT_BOUNDS
    lon_vertices.attrs = cf.coords[cf.LON_BOUNDS]
    lat_vertices.attrs = cf.coords[cf.LAT_BOUNDS]
    return xr.merge([lat_vertices, lon_vertices])
