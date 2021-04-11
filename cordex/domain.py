# -*- coding: utf-8 -*-
# flake8: noqa
"""Domain module

This module defines preconfigured CORDEX domain. The :class:`Domain` class
is a kind of wrapper for the :class:`Grid` class to connect a grid
with meta information and easy to use functions that work on the member grid.

Domains can either be defined statically in this module or read from a csv table
that should be defined in the :mod:`table`.

Example:

    To get a list of available implementations, create cordex domains, write
    them to netcdf with some dummy data, you can use ,e.g.,::

        from cordex import domain as dm

        for short_name in dm.domains():
            print('creating domain: {}'.format(short_name))
            domain = dm.domain(short_name)
            domain.to_netcdf(short_name+'.nc', dummy='topo')

"""

import numpy as np
import pandas as pd
import xarray as xr

from .tables import domains as TABLES
from . import tables
from . import cf
from . import utils


def domain_from_table(short_name, table):
    """creates domain instance from a pandas dataframe row.
    """
    return Domain(short_name=short_name, **dict(table.loc[short_name]))


def domain(short_name, dummy=False, **kwargs):
    """Creates an xarray dataset containg the domain grid definitions.

    Parameters
    ----------
    short_name:
        Name of the Cordex Domain.
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be 
        used to create some dummy topography data. dummy data is useful for 
        looking at the domain with ncview.

    Returns
    -------
    Dataset : xarray.core.Dataset
        Dataset containing the coordinates.
        
    """
    config = pd.concat(TABLES.values()).loc[short_name]
    print(config)
    return create_dataset(**config, dummy=dummy)


def create_dataset(nlon, nlat, dlon, dlat, ll_lon, ll_lat, 
        pollon, pollat, dummy=False,  **kwargs):
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


def _get_dataset(rlon, rlat, lon, lat, pole, dummy=None, mapping_name=None, attrs=None):
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    data_vars={mapping_name: pole}

    ds = xr.Dataset(
        data_vars=data_vars,
        coords=dict(
            rlon=(["rlon"], rlon),
            rlat=(["rlat"], rlat),
            lon=(["rlat", "rlon"], lon),
            lat=(["rlat", "rlon"], lat),
        ),
        attrs=attrs
    )

    for key, coord  in ds.coords.items():
        coord.encoding['_FillValue'] = None 
        coord.attrs = cf.coords[key]

    if dummy:
        if dummy is True:
            dummy_name = 'dummy'
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
        dummy.attrs = { 'grid_mapping' : mapping_name,
                        'coordinates'  :  'lon lat' }
        ds[dummy_name] = dummy
        if dummy_name == 'topo':
            from cdo import Cdo
            tmp = utils.get_tempfile()
            ds.to_netcdf(tmp)
            topo = Cdo().topo(tmp, returnXDataset=True)['topo'][:]
            ds[dummy_name] = topo

    return ds


def _grid_mapping(pollon, pollat, mapping_name=None):
    """creates a grid mapping DataArray object"""
    if mapping_name is None:
        mapping_name = cf.DEFAULT_MAPPING_NCVAR
    da = xr.DataArray(np.zeros((), dtype=np.int32))
    attrs = cf.mapping.copy()
    attrs['grid_north_pole_longitude'] = pollon 
    attrs['grid_north_pole_latitude']  = pollat 
    da.attrs = attrs
    return da


def _init_grid(nlon, nlat, dlon, dlat, ll_lon, ll_lat):
    """create coordinate arrays from lower left longitude and latitude"""
    rlon = np.array([ll_lon+i*dlon for i in range(0,nlon)], dtype=np.float64)
    rlat = np.array([ll_lat+i*dlat for i in range(0,nlat)], dtype=np.float64)
    return rlon, rlat


def _stack(x, y):
    """Stack 1d arrays into 2d fields."""
    tmp_x = np.array(x).squeeze()
    tmp_y = np.array(y).squeeze()
    return np.vstack(len(tmp_y)*(tmp_x,)), np.hstack(len(tmp_x)*(tmp_y[:, np.newaxis],))


def rotated_coord_transform(lon, lat, np_lon, np_lat,
                            direction='rot2geo'):
    """Transforms a coordinate into a rotated grid coordinate and vice versa.

    The coordinates have to given in degree and will be returned in degree.

    Parameters
    ----------
    lon : float
        Longitude coordinate.
    lat : float
        Latitude coordinate.
    np_lon : float
        Longitude coordinate of the rotated pole.
    np_lat : float
        Latitude coordinate of the rotated pole.
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
    lon = (lon * np.pi) / 180.
    lat = (lat * np.pi) / 180.

    theta = 90. - np_lat # Rotation around y-axis
    phi = np_lon + 180.  # Rotation around z-axis

    # Convert degrees to radians
    phi = (phi * np.pi) / 180.
    theta = (theta * np.pi) / 180.

    # Convert from spherical to cartesian coordinates
    x = np.cos(lon) * np.cos(lat)
    y = np.sin(lon) * np.cos(lat)
    z = np.sin(lat)

    # Regular -> Rotated
    if direction == 'geo2rot':

        x_new = (np.cos(theta) * np.cos(phi) * x +
                 np.cos(theta) * np.sin(phi) * y +
                 np.sin(theta) * z)
        y_new = (- np.sin(phi) * x +
                   np.cos(phi) * y)
        z_new = (- np.sin(theta) * np.cos(phi) * x -
                   np.sin(theta) * np.sin(phi) * y +
                   np.cos(theta) * z)

    # Rotated -> Regular
    elif direction == 'rot2geo':
    
        phi = - phi
        theta = - theta
    
        x_new = (np.cos(theta) * np.cos(phi) * x +
                 np.sin(phi) * y +
                 np.sin(theta) * np.cos(phi) * z)
        y_new = (- np.cos(theta) * np.sin(phi) * x +
                   np.cos(phi) * y -
                   np.sin(theta) * np.sin(phi) * z)
        z_new = (- np.sin(theta) * x +
                   np.cos(theta) * z)

    # Convert cartesian back to spherical coordinates
    lon_new = np.arctan2(y_new, x_new)
    lat_new = np.arcsin(z_new)

    # Convert radians back to degrees
    lon_new = (lon_new * 180.) / np.pi
    lat_new = (lat_new * 180.) / np.pi;

    return lon_new, lat_new


