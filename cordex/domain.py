from warnings import warn

import cf_xarray as cfxr  # noqa
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import CRS

from . import cf
from .config import nround
from .tables import domains
from .transform import grid_mapping, transform, transform_bounds
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
    bounds=False,
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

        .. deprecated:: v0.5.0
            The ``add_vertices`` parameter is deprecated in favor
            of the ``bounds`` parameter, and will be removed in a future
            version.

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
    bounds: bool
        Add spatial bounds to longitude and latitude coordinates.

        .. versionadded:: v0.5.0

    Returns
    -------
    Grid : xr.Dataset
        Dataset containing a CORDEX grid.

    References
    ----------
    Please refer to the CF conventions document : https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#grid-mappings-and-projections

    Example
    -------

    To create a cordex rotated pole domain dataset, you can use ,e.g.,::

        import cordex as cx

        eur11 = cx.cordex_domain('EUR-11')

    """
    if add_vertices is True:
        warn(
            "add_vertices keyword is deprecated, please use the bounds keyword instead",
            DeprecationWarning,
            stacklevel=2,
        )
        bounds = True

    if attrs is None:
        attrs = {}
    if tables is None:
        tables = domains.table
    if isinstance(tables, list):
        tables = pd.concat(tables)
    config = tables.replace(np.nan, None).loc[short_name]

    return create_dataset(
        **config,
        short_name=short_name,
        dummy=dummy,
        add_vertices=False,
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
    short_name=None,
    dummy=False,
    add_vertices=False,
    attrs=None,
    mapping_name=None,
    bounds=False,
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
    short_name : str
        CORDEX domain identifier, goes into the ``CORDEX_domain`` global attribute.
    dummy : str or logical
        Name of dummy field, if dummy=topo, the cdo topo operator will be
        used to create some dummy topography data. dummy data is useful for
        looking at the domain with ncview.
    add_vertices : bool
        Add grid boundaries in the global coordinates (lon_vertices and lat_vertices).

        .. deprecated:: v0.5.0
            The ``add_vertices`` parameter is deprecated in favor
            of the ``bounds`` parameter, and will be removed in a future
            version.

    attrs: str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name: str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    bounds: bool
        Add spatial bounds to longitude and latitude coordinates.

        .. versionadded:: v0.5.0

    Returns
    -------
    Grid : Dataset
        Dataset containing a CORDEX grid.

    """
    if add_vertices is True:
        warn(
            "add_vertices keyword is deprecated, please use the bounds keyword instead",
            DeprecationWarning,
            stacklevel=2,
        )
        bounds = True

    rotated = True
    if attrs == "CORDEX":
        attrs = cf.DEFAULT_CORDEX_ATTRS
    elif attrs is None:
        attrs = {}
    if name:
        attrs["CORDEX_domain"] = name
    # remove inconsistencies in keyword names
    if short_name:
        attrs["CORDEX_domain"] = short_name
    if pollon is None or pollat is None:
        rotated = False
    try:
        if np.isnan(pollon) or np.isnan(pollat):
            rotated = False
    except Exception:
        pass

    x, y = _lin_coord(nlon, dlon, ll_lon), _lin_coord(nlat, dlat, ll_lat)

    if rotated is True:
        ds = _get_rotated_dataset(
            x,
            y,
            #  lon,
            #  lat,
            pollon,
            pollat,
            bounds=bounds,
            mapping_name=mapping_name,
            attrs=attrs,
        )
    else:
        ds = _get_regular_dataset(
            x,
            y,
            bounds=bounds,
            attrs=attrs,
        )

    if dummy:
        ds = _add_dummy(ds, dummy)

    return ds


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
    elif isinstance(tables, list):
        tables = pd.concat(tables)

    config = tables.replace(np.nan, None).loc[short_name]
    # return config
    return {**{"short_name": short_name}, **config.to_dict()}


def _get_regular_dataset(
    lon,
    lat,
    bounds=False,
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

    if bounds is True:
        ds = ds.cf.add_bounds(("lon", "lat"))

    return ds


def _get_rotated_dataset(
    rlon,
    rlat,
    pollon,
    pollat,
    bounds=False,
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

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cf.coords[key]

    if bounds is True:
        ds = transform_bounds(ds)
        # ds = xr.merge([ds, v])
        ds[cf.LON_NAME].attrs["bounds"] = cf.LON_BOUNDS
        ds[cf.LAT_NAME].attrs["bounds"] = cf.LAT_BOUNDS

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
    warn(
        "vertices is deprecated, please use transform_bounds instead",
        DeprecationWarning,
        stacklevel=2,
    )
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


def _crop_to_domain(ds, domain_id, drop=True):
    domain = cordex_domain(domain_id)
    x_mask = ds.rlon.round(8).isin(domain.cf["X"])
    y_mask = ds.rlat.round(8).isin(domain.cf["Y"])
    return ds.where(x_mask & y_mask, drop=drop)
