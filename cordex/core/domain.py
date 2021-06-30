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

from ..domains import read_cordex_tables
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
        return list(read_cordex_tables()[table_name].index)
    else:
        return list(pd.concat(read_cordex_tables().values()).index)


def cordex_domain(short_name, dummy=False, tables=None):
    """Creates an xarray dataset containg the domain grid definitions.

    Parameters
    ----------
    short_name:
        Name of the Cordex Domain.
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.
    tables: dataframe or list of dataframes, default: cordex_tables
        Tables from which to look up the grid information. Index in the table
        should be the short name of the domain, e.g., `EUR-11`. If no table is
        provided, all standard tables will be searched.

    Returns
    -------
    Dataset : xarray.core.Dataset
        Dataset containing the coordinates.

    """
    if tables is None:
        tables = list(read_cordex_tables().values())
    if len(tables) > 1:
        config = pd.concat(tables).loc[short_name]
    else:
        config = tables.loc[short_name]
    return create_dataset(**config, dummy=dummy)


def create_dataset(
    nlon, nlat, dlon, dlat, ll_lon, ll_lat, pollon, pollat, dummy=False, **kwargs
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
    """
    rlon, rlat = _init_grid(nlon, nlat, dlon, dlat, ll_lon, ll_lat)
    lon, lat = rotated_coord_transform(*_stack(rlon, rlat), pollon, pollat)
    pole = _grid_mapping(pollon, pollat)
    return _get_dataset(rlon, rlat, lon, lat, pole, dummy=dummy)


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
        tables = list(read_cordex_tables().values())
    if len(tables) > 1:
        config = pd.concat(tables).loc[short_name]
    else:
        config = tables.loc[short_name]
    return {**{"short_name": short_name}, **dict(**config)}


def _get_dataset(rlon, rlat, lon, lat, pole, dummy=None, mapping_name=None, attrs=None):
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    data_vars = {mapping_name: pole}

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            rlon=(["rlon"], rlon),
            rlat=(["rlat"], rlat),
            lon=(["rlat", "rlon"], lon),
            lat=(["rlat", "rlon"], lat),
        ),
        attrs=attrs,
    )

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cf.coords[key]

    if dummy:
        if dummy is True:
            dummy_name = "dummy"
        else:
            dummy_name = dummy
        dummy = xr.DataArray(
            data=np.zeros((len(rlat), len(rlon))),
            dims=["rlat", "rlon"],
            coords=dict(
                rlon=(["rlon"], rlon),
                rlat=(["rlat"], rlat),
                lon=(["rlat", "rlon"], lon),
                lat=(["rlat", "rlon"], lat),
            ),
        )
        dummy.attrs = {"grid_mapping": mapping_name, "coordinates": "lon lat"}
        ds[dummy_name] = dummy
        if dummy_name == "topo":
            from cdo import Cdo

            tmp = utils.get_tempfile()
            ds.to_netcdf(tmp)
            topo = Cdo().topo(tmp, returnXDataset=True)["topo"][:]
            ds[dummy_name] = topo

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
