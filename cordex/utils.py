import tempfile

import numpy as np
import pandas as pd
from pyproj import Transformer, Proj
from .tables import domains


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


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


def _cell_area(ds, R=6371000):
    """Compute cell area of a regular spherical grid.

    Parameters
    ----------
    ds : str
        Dataset containing longitude and latitude coordinates.
    R : float
        Earth radius in units [m]. Defaults to 6371000 meters.
    """

    dphi, dtheta = (
        np.deg2rad(
            ds.cf[c]
            # pad values to account for differentiation
            .pad({ds.cf[c].dims[0]: (1, 0)}, mode="reflect", reflect_type="odd").diff(
                ds.cf[c].dims[0]
            )
        )
        for c in ("X", "Y")
    )
    dOmega = np.cos(np.deg2rad(ds.cf["Y"])) * dtheta * dphi
    return R**2 * dOmega


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


def cell_area(ds, R=6371000, attrs=True):
    r"""Compute cell areas for a regular spherical grid.

    Parameters
    ----------
    ds : str
        Dataset containing a regular grid with longitude and latitude
        coordinates that can be understood by cf_xarray.
    R : float
        Earth radius in units [m]. Defaults to 6371000 meters.
    attrs: logical or str
        If True, add attributes for grid-cell area. If ``"CF"``,
        add CF attributes for atmospheric grid-cell area.

    Returns
    -------
    Cell area : xr.DataArray
        DataArray containg the size of each grid cell in units [m2]

    Notes
    -----
    The solid angle differential of the sphere is computed as

    .. math::
        d\Omega = \sin\theta\,d\theta\,d\phi

    with the surface element

    .. math::
        dA = R^2 d\Omega

    References
    ----------

    https://en.wikipedia.org/wiki/Solid_angle

    """

    da = _cell_area(ds, R)

    if attrs:
        da.name = "areacell"
        da.attrs = {
            "standard_name": "cell_area",
            "units": "m2",
            "long_name": "Grid-Cell Area",
        }
    if attrs == "CF":
        da.name = "areacella"
        da.attrs.update(
            {
                "cell_methods": "area: sum",
                "cell_measures": "area: areacella",
                "long_name": "Atmosphere Grid-Cell Area",
            }
        )
    if not attrs:
        da.attrs = {}

    return da


def create_polygon(obj, crs=None):
    """create polygon in rotated pole coords"""
    from shapely.geometry import Polygon

    coords = [
        [obj.cf["X"].min(), obj.cf["Y"].min()],
        [obj.cf["X"].max(), obj.cf["Y"].min()],
        [obj.cf["X"].max(), obj.cf["Y"].max()],
        [obj.cf["X"].min(), obj.cf["Y"].max()],
    ]
    return Polygon(coords)


def transform_polygon(polygon, crs, segmentize=None):
    """
    Segmentize and transform a polygon to latitude/longitude coordinates.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        The polygon to transform.
    crs : pyproj.CRS or str
        The source coordinate reference system.
    segmentize : float, optional
        The maximum segment length for segmentizing the polygon.

    Returns
    -------
    shapely.geometry.Polygon
        The transformed polygon in latitude/longitude coordinates.
    """
    from shapely.ops import transform

    transformer = Transformer.from_proj(crs, Proj(init="epsg:4326"))
    if segmentize:
        polygon = polygon.segmentize(segmentize)
    return transform(transformer.transform, polygon)


def get_grid_mapping_attrs(domain_id):
    """
    Get grid mapping attributes for a given domain.

    Parameters
    ----------
    domain_id : str
        The domain identifier.

    Returns
    -------
    dict
        A dictionary containing grid mapping attributes.
    """
    grid = domain_info(domain_id)
    attrs = grid.copy()
    attrs["grid_north_pole_longitude"], attrs["grid_north_pole_latitude"] = (
        grid["pollon"],
        grid["pollat"],
    )
    return attrs
