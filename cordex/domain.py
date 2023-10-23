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


def _locate_domain_id(domain_id, table):
    """Locate domain_id in domain table trying different indexes.

    First, it is assumed that domain_id can be found in the tables index.
    If it is not found in the index, a number of different colums are
    tried as index (``short_name``, ``domain_id``, ``CORDEX_domain``).

    """

    indexes = [table.index.name]
    # additional indexes to try
    indexes.extend(["short_name", "domain_id", "CORDEX_domain"])
    # removed duplicates
    indexes = list(dict.fromkeys(indexes))

    table = table.replace(np.nan, None)

    for i in indexes:
        if i in table.reset_index().columns:
            if domain_id in table.reset_index()[i].values:
                return (
                    table.reset_index()
                    .set_index(i)
                    .loc[[domain_id]]
                    .reset_index()
                    .iloc[0]
                )

    return table.replace(np.nan, None).loc[domain_id]


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
    domain_id,
    dummy=False,
    add_vertices=False,
    tables=None,
    attrs=None,
    mapping_name=None,
    bounds=False,
    mip_era="CMIP5",
):
    """Creates an xarray dataset containg the domain grid definitions.

    Parameters
    ----------
    domain_id : str
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

    tables : dataframe or list of dataframes, default: cordex_tables
        Tables from which to look up the grid information. Index in the table
        should be the short name of the domain, e.g., `EUR-11`. If no table is
        provided, all standard tables will be searched.
    attrs : str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name : str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    bounds : bool
        Add spatial bounds to longitude and latitude coordinates.

        .. versionadded:: v0.5.0
    mip_era : str
        The mip_era keyword determines the vocabulary for dimensions, coordinates and
        attributes.

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

    config = _locate_domain_id(domain_id, tables)

    return create_dataset(
        **config,
        # domain_id=domain_id,
        dummy=dummy,
        add_vertices=False,
        attrs=attrs,
        mapping_name=mapping_name,
        bounds=bounds,
        mip_era=mip_era,
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
    domain_id=None,
    dummy=False,
    add_vertices=False,
    attrs=None,
    mapping_name=None,
    bounds=False,
    mip_era="CMIP5",
    **kwargs,
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
    domain_id : str
        Domain identifier, goes into the ``CORDEX_domain`` or ``domain_id`` global attribute.
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

    attrs : str or dict
        Global attributes that should be added to the dataset. If `attrs='CORDEX'`
        a set of standard CF global attributes.
    mapping_name : str
        Variable name of the grid mapping, if mapping_name is `None`, the CF standard
        variable name is used.
    bounds : bool
        Add spatial bounds to longitude and latitude coordinates.

        .. versionadded:: v0.5.0
    mip_era : str
        The mip_era keyword determines the vocabulary for dimensions, coordinates and
        attributes.

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

    cv = cf.vocabulary[mip_era]

    rotated = True
    if attrs == "CORDEX":
        attrs = cv["default_global_attrs"]
    elif attrs is None:
        attrs = {}

    if name:
        attrs[cv["domain_id"]] = name
    # remove inconsistencies in keyword names
    if domain_id:
        attrs[cv["domain_id"]] = domain_id
    if cv["domain_id"] in kwargs:
        attrs[cv["domain_id"]] = kwargs[cv["domain_id"]]

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
            cv=cv,
        )
    else:
        ds = _get_regular_dataset(
            x,
            y,
            bounds=bounds,
            attrs=attrs,
            cv=cv,
        )

    if dummy:
        ds = _add_dummy(ds, dummy)

    return ds


def domain_info(domain_id, tables=None):
    """Returns a dictionary containg the domain grid definitions.

    Returns a dictionary with grid information according to the
    Cordex archive specifications.

    See https://is-enes-data.github.io/cordex_archive_specifications.pdf

    Parameters
    ----------
    domain_id:
        Cordex domain identifier.

    Returns
    -------
    domain info : dict
        Dictionary containing the grid information.

    """
    if tables is None:
        tables = domains.table
    elif isinstance(tables, list):
        tables = pd.concat(tables)

    return _locate_domain_id(domain_id, tables).to_dict()


def _get_regular_dataset(
    x,
    y,
    bounds=False,
    attrs=None,
    cv=None,
):
    xdim = cv["dims"]["LON"]
    ydim = cv["dims"]["LAT"]
    ds = xr.Dataset(
        data_vars=None,
        coords={
            xdim: (xdim, x),
            ydim: (ydim, y),
        },
        attrs=attrs,
    )

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cv["coords"][key]

    ds[xdim].attrs["axis"] = "X"
    ds[ydim].attrs["axis"] = "Y"

    if bounds is True:
        ds = ds.cf.add_bounds((xdim, ydim))

    return ds


def _get_rotated_dataset(
    x,
    y,
    pollon,
    pollat,
    bounds=False,
    mapping_name=None,
    attrs=None,
    cv=None,
):
    xdim = cv["dims"]["X"]
    ydim = cv["dims"]["Y"]
    lon_dim = cv["dims"]["LON"]
    lat_dim = cv["dims"]["LAT"]
    lon_bounds = cv["dims"]["LON_BOUNDS"]
    lat_bounds = cv["dims"]["LAT_BOUNDS"]
    mapping_name = mapping_name or cv["default_mapping_ncvar"]

    mapping = grid_mapping(pollon, pollat, mapping_name)
    # lon, lat = transform(x, y, src_crs=CRS.from_cf(mapping.attrs))

    ds = xr.Dataset(
        data_vars={mapping.name: mapping},
        coords={
            xdim: (xdim, x),
            ydim: (ydim, y),
        },
        attrs=attrs,
    )

    lon, lat = transform(ds[xdim], ds[ydim], src_crs=CRS.from_cf(mapping.attrs))
    ds = ds.assign_coords({lon_dim: lon, lat_dim: lat})

    for key, coord in ds.coords.items():
        coord.encoding["_FillValue"] = None
        coord.attrs = cv["coords"][key]

    if bounds is True:
        ds = transform_bounds(ds, trg_dims=(lon_bounds, lat_bounds))
        # ds = xr.merge([ds, v])
        ds[lon_dim].attrs["bounds"] = lon_bounds
        ds[lat_dim].attrs["bounds"] = lat_bounds

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
    x_mask = ds.cf["X"].round(8).isin(domain.cf["X"])
    y_mask = ds.cf["Y"].round(8).isin(domain.cf["Y"])
    return ds.where(x_mask & y_mask, drop=drop)
