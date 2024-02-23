import tempfile

import numpy as np


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


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
