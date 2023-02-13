"""CORDEX Cmorization utilities.
"""
import datetime as dt
import json
from warnings import warn

import cftime as cfdt
import numpy as np
import xarray as xr

from .. import cordex_domain

xr.set_options(keep_attrs=True)

loffsets = {"3H": dt.timedelta(hours=1, minutes=30), "6H": dt.timedelta(hours=3)}


def _get_loffset(time):
    return loffsets.get(time, None)


# def ensure_cftime(func):
#    def wrapper(date, **kwargs):
#        return func(_to_cftime(date), **kwargs)
#
#    return wrapper


# def to_cftime(date, calendar="gregorian"):
#     """Convert datetime object to cftime object.

#     Parameters
#     ----------
#     date : datetime object
#         Datetime object.
#     calendar : str
#         Calendar of the cftime object.

#     Returns
#     -------
#     cftime : cftime object
#         Cftime ojbect.

#     """
#     if type(date) == dt.date:
#         date = dt.datetime.combine(date, dt.time())
#     elif isinstance(date, cfdt.datetime):
#         # do nothing
#         return date
#     return cfdt.datetime(
#         date.year,
#         date.month,
#         date.day,
#         date.hour,
#         date.minute,
#         date.second,
#         date.microsecond,
#         calendar=calendar,
#     )


def to_cftime(date, calendar="standard"):
    """Convert date to cftime object

    Can handle all CMIP6 calendars.

    Parameters
    ----------
    date : datetime object, str
        Input date.
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
    elif isinstance(date, str):
        # xarray hack for cftime.strptime
        return xr.cftime_range(start=date, end=date, calendar=calendar)[0]
        # date = pd.to_datetime(date)
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
    except Exception:
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


def _month_bounds(date):
    """Determine the bounds of the current month.

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
    begin = date.replace(day=1, hour=0, minute=0, second=0)
    # this does not work with cftime
    # end = (date + reld.relativedelta(months=1)).replace(day=1)
    if month == 12:
        year = date.year + 1
        month = 1
    else:
        year = date.year
        month = date.month + 1
    end = date.replace(day=1, year=year, month=month, hour=0, minute=0, second=0)
    return begin, end


def _mid_of_month(date, return_bounds=False):
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
    bounds = _month_bounds(date)
    mid = bounds[0] + 0.5 * (bounds[1] - bounds[0])
    if return_bounds is True:
        return bounds
    return mid


def _helper(date):
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
    return np.array([_month_bounds(d) for d in date], dtype="o,o")


def month_bounds(ds):
    return xr.apply_ufunc(
        _helper,
        ds.time,
        # input_core_dims=[["time"]],
        # output_core_dims=[[],[]],
        # input_core_dims=[[]],
        # output_core_dims=[[]],
        vectorize=False,
    )


def mid_of_month(ds, add_bounds=False):
    """Determine the mid of the current month.

    Parameters
    ----------
    ds : Dataset or DataArray
        Date in the current month.

    Returns
    -------
    mid_of_month : datetime object
        Mid date of the current month.

    """
    ds = ds.copy(deep=False)
    if add_bounds is True:
        return xr.apply_ufunc(
            _mid_of_month,
            ds.time,
            True,
            input_core_dims=[[], []],
            output_core_dims=[[], []],
            vectorize=True,
        )
        # output_core_dims=[[ds.time.dims],[ds.time.dims]],vectorize=True)
    return xr.apply_ufunc(_mid_of_month, ds.time, False, vectorize=True)


def _get_pole(ds):
    """returns the first pole we find in the dataset"""
    pol_names = ["rotated_latitude_longitude", "rotated_pole"]
    for pol in pol_names:
        if pol in ds:
            return ds[pol]
    warn("no grid_mapping found in dataset, tried: {}".format(pol_names))
    return None


def _get_grid_definitions(CORDEX_domain, **kwargs):
    return cordex_domain(CORDEX_domain, add_vertices=True, **kwargs)


def _get_cordex_pole(CORDEX_domain):
    return cordex_domain(CORDEX_domain).rotated_latitude_longitude


def _encode_time(time):
    """encode xarray time axis into cf values

    see https://github.com/pydata/xarray/issues/4412

    """
    return xr.conventions.encode_cf_variable(time)


def _read_cmor_table(table):
    return _read_json_file(table)


def _read_json_file(filename):
    with open(filename) as f:
        data = json.load(f)
    return data


def _get_cfvarinfo(cf_varname, table):
    data = _read_cmor_table(table)
    return data["variable_entry"].get(cf_varname, None)


def _get_time_cell_method(cf_varname, table):
    return _strip_time_cell_method(_get_cfvarinfo(cf_varname, table))


def _strip_time_cell_method(cfvarinfo):
    try:
        return cfvarinfo["cell_methods"].split("time:")[1].strip()
    except Exception:
        return None
