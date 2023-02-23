import tempfile

import numpy as np

import cordex as cx

from .config import nround


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


def _get_info(ds, tables=None, precision=nround):
    if tables is None:
        tables = cx.domains.table.replace(np.nan, None)
    try:
        x = ds.cf["X"]
        y = ds.cf["Y"]
    except KeyError:
        x = ds.rlon
        y = ds.rlat
    nlon = x.size
    nlat = y.size
    dlon = x.diff(x.dims[0]).values.item(0)
    dlat = y.diff(y.dims[0]).values.item(0)
    ll_lon = x[0].min().values.item(0)
    ll_lat = y[0].min().values.item(0)
    ur_lon = x[-1].min().values.item(0)
    ur_lat = y[-1].min().values.item(0)
    try:
        pollon = ds.cf["grid_mapping"].grid_north_pole_longitude
        pollat = ds.cf["grid_mapping"].grid_north_pole_latitude
    except KeyError:
        pollon = None
        pollat = None
    coords = {
        "nlon": nlon,
        "nlat": nlat,
        "ll_lon": ll_lon,
        "ur_lon": ur_lon,
        "ll_lat": ll_lat,
        "ur_lat": ur_lat,
        "dlon": dlon,
        "dlat": dlat,
        "pollon": pollon,
        "pollat": pollat,
    }
    # round
    info = {
        k: (np.round(v, nround) if isinstance(v, float) else v)
        for k, v in coords.items()
    }
    return info


def _guess_domain(ds, tables=None):
    if tables is None:
        tables = cx.domains.table  # .replace(np.nan, None)
    try:
        info = _get_info(ds, tables)
    except Exception as e:
        print(e)
        raise Exception(
            "Could not determine domain, only rotated_latitude_longitude supported."
        )
    filt = tables
    for k, v in info.items():
        if filt.empty:
            return None
        if v:
            filt = filt[np.isclose(filt[k], v)]
        else:
            filt = filt[np.isnan(filt[k])]  # | filt[k] is None]
    # reset index and convert to dict
    return filt.reset_index().iloc[0].replace(np.nan, None).to_dict()
