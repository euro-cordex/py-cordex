"""CORDEX Cmorization utilities.
"""
import numpy as np
import xarray as xr
import datetime as dt
import cftime as cfdt
from dateutil import relativedelta as reld

from ..core.domain import map_crs

xr.set_options(keep_attrs=True)

loffsets = {"3H": dt.timedelta(hours=1, minutes=30), "6H": dt.timedelta(hours=3)}


def _get_loffset(time):
    return loffsets.get(time, None)


def ensure_cftime(func):
    def wrapper(date, **kwargs):
        return func(_to_cftime(date), **kwargs)

    return wrapper


def to_cftime(date, calendar="proleptic_gregorian"):
    """Convert datetime object to cftime object.

    Parameters
    ----------
    date : datetime object
        Datetime object.
    calendar : str
        Calendar of the cftime object.

    Returns
    -------
    cftime : cftime object
        Cftime ojbect.

    """
    if type(date) == dt.date:
        date = dt.datetime.combine(date, dt.time())
    elif isinstance(date, cfdt.datetime):
        # do nothing
        return date
    return cfdt.datetime(
        date.year,
        date.month,
        date.day,
        date.hour,
        date.minute,
        date.second,
        date.microsecond,
        calendar=calendar,
    )


def _seasons_bounds(year, calendar=None):
    if calendar is None:
        import datetime as dt

        args = {}
    else:
        # calendar requires cftime
        import cftime as dt

        args = {"calendar": calendar}
    return {
        "DJF": (dt.datetime(year - 1, 12, 1, **args), dt.datetime(year, 3, 1, **args)),
        "MAM": (dt.datetime(year, 3, 1, **args), dt.datetime(year, 6, 1, **args)),
        "JJA": (dt.datetime(year, 6, 1, **args), dt.datetime(year, 9, 1, **args)),
        "SON": (dt.datetime(year, 9, 1, **args), dt.datetime(year, 12, 1, **args)),
    }


def season_bounds(date):
    """Determines the temporal bounds of the meteorological season.

    Uses the month to determine the season and returns the
    temporal bounds of the season.

    Parameters
    ----------
    date : datetime object
        Date in the current season.

    Returns
    -------
    season : tuple of datetime objects
        Temporal bounds of the current meteorological season.

    """
    month = date.month
    if month != 12:
        year = date.year
    else:
        year = date.year + 1
    try:
        calendar = date.calendar
    except:
        calendar = None
    seasons_bounds = _seasons_bounds(year, calendar=calendar)
    return seasons_bounds[season(date)]


def _seasons():
    seasons = [
        ("DJF", (12, 1, 2)),
        ("MAM", (3, 4, 5)),
        ("JJA", (6, 7, 8)),
        ("SON", (9, 10, 11)),
    ]
    return seasons


# @ensure_cftime
def season(date):
    """Determines the meteorological season.

    Uses the month to determine the season.

    Parameters
    ----------
    date : datetime object
        Date in the current season.

    Returns
    -------
    season : str
        Meteorological season of the current date.

    """
    return next(season for season, months in _seasons() if date.month in months)


def mid_of_season(date):
    """Determine the mid of the current season

    Parameters
    ----------
    date : datetime object
        Date in the current season.

    Returns
    -------
    mid_of_season : datetime object
        Mid date of the current season.

    """
    bounds = season_bounds(date)
    return bounds[0] + 0.5 * (bounds[1] - bounds[0])


def month_bounds(date):
    """Determine the mid of the current month.

    Parameters
    ----------
    date : datetime object
        Date in the current month.

    Returns
    -------
    month_bounds : tuple of datetime object
        Temporal bounds of the current month.

    """
    if type(date) == dt.date:
        date = dt.datetime.combine(date, dt.time())
    month = date.month
    begin = date.replace(day=1)
    # this does not work with cftime
    # end = (date + reld.relativedelta(months=1)).replace(day=1)
    if month == 12:
        year = date.year + 1
        month = 1
    else:
        year = date.year
        month = date.month + 1
    end = date.replace(day=1, year=year, month=month)
    return (begin, end)


def mid_of_month(date):
    """Determine the mid of the current month.

    Parameters
    ----------
    date : datetime object
        Date in the current month.

    Returns
    -------
    mid_of_month : datetime object
        Mid date of the current month.

    """
    bounds = month_bounds(date)
    return bounds[0] + 0.5 * (bounds[1] - bounds[0])


def _dcoord(coord):
    dcoord = coord.values[1:] - coord.values[:-1]
    dcoord = np.insert(dcoord, 0, dcoord[0])
    return dcoord


def _bounds(coord):
    dc = _dcoord(coord)
    left = coord - 0.5 * dc
    right = coord + 0.5 * dc
    left.name = "left"
    right.name = "right"
    return xr.merge([left, right])


def vertices(rlon, rlat, src_crs, trg_crs=None):
    """Compute lon and lat vertices.

    Transformation of rlon vertices and rlat vertices
    into the target crs.

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
    lon_vertices = xr.concat([v1[0], v2[0], v3[0], v4[0]], dim="vertices").transpose(
        ..., "vertices"
    )
    lat_vertices = xr.concat([v1[1], v2[1], v3[1], v4[1]], dim="vertices").transpose(
        ..., "vertices"
    )
    lon_vertices.name = "lon_vertices"
    lat_vertices.name = "lat_vertices"
    return xr.merge([lon_vertices, lat_vertices])
