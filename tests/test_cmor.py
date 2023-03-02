import datetime as dt

import cftime
import cftime as cfdt
import numpy as np
import pytest
import xarray as xr

import cordex as cx
from cordex import cmor
from cordex.cmor import utils


def test_cfdt():
    assert cmor.to_cftime(dt.datetime(2000, 1, 1, 1)) == cfdt.datetime(2000, 1, 1, 1)
    assert cmor.to_cftime(dt.date(2000, 1, 1)) == cfdt.datetime(2000, 1, 1)
    assert cmor.to_cftime("2000-01-01T01:00:00") == cfdt.datetime(2000, 1, 1, 1)
    assert cmor.to_cftime("2000-02-30T00:00:00", calendar="360_day") == cfdt.datetime(
        2000, 2, 30, calendar="360_day"
    )


@pytest.mark.parametrize("dt", [dt, cfdt])
def test_season(dt):
    assert cmor.season(dt.datetime(2000, 1, 1)) == "DJF"
    assert cmor.season(dt.datetime(2001, 3, 1)) == "MAM"
    assert cmor.season(dt.datetime(2001, 6, 1)) == "JJA"
    assert cmor.season(dt.datetime(2001, 9, 1)) == "SON"
    assert cmor.season(dt.datetime(2001, 12, 1)) == "DJF"
    bounds = (dt.datetime(1999, 12, 1, 0, 0), dt.datetime(2000, 3, 1, 0, 0))
    assert cmor.season_bounds(dt.datetime(2000, 1, 1)) == bounds
    assert cmor.mid_of_season(dt.datetime(2000, 1, 1)) == dt.datetime(
        2000, 1, 15, 12, 0
    )


def test_month():
    assert cmor.month_bounds(dt.datetime(2000, 1, 10)) == (
        dt.datetime(2000, 1, 1),
        dt.datetime(2000, 2, 1),
    )
    assert cmor.month_bounds(dt.datetime(2000, 12, 10)) == (
        dt.datetime(2000, 12, 1),
        dt.datetime(2001, 1, 1),
    )
    assert cmor.mid_of_month(dt.datetime(2000, 1, 1)) == dt.datetime(2000, 1, 16, 12)
    # leap years
    assert cmor.mid_of_month(dt.datetime(2000, 2, 1)) == dt.datetime(2000, 2, 15, 12)
    assert cmor.mid_of_month(dt.datetime(2001, 2, 1)) == dt.datetime(2001, 2, 15)


def test_cfmonth():
    assert cmor.month_bounds(dt.datetime(2000, 1, 10)) == (
        dt.datetime(2000, 1, 1),
        dt.datetime(2000, 2, 1),
    )
    # leap years
    assert cmor.mid_of_month(cfdt.datetime(2000, 2, 1)) == cfdt.datetime(
        2000, 2, 15, 12
    )
    assert cmor.mid_of_month(cfdt.datetime(2001, 2, 1)) == cfdt.datetime(2001, 2, 15)
    assert cmor.mid_of_month(
        cfdt.datetime(2000, 2, 1, calendar="360_day")
    ) == cfdt.datetime(2000, 2, 16, calendar="360_day")
    assert cmor.mid_of_month(
        cfdt.datetime(2001, 2, 1, calendar="360_day")
    ) == cfdt.datetime(2001, 2, 16, calendar="360_day")


def test_mid_of_month():
    time_axis = xr.DataArray(
        xr.cftime_range("2005-01", periods=12, freq="MS"), dims="time"
    )

    expect = np.array(
        [
            cfdt.DatetimeGregorian(2005, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 2, 15, 0, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 3, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 4, 16, 0, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 5, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 6, 16, 0, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 7, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 8, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 9, 16, 0, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 10, 16, 12, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 11, 16, 0, 0, 0, 0, has_year_zero=False),
            cfdt.DatetimeGregorian(2005, 12, 16, 12, 0, 0, 0, has_year_zero=False),
        ]
    )

    expect = xr.DataArray(expect, dims="time")

    mid = utils.mid_of_month(time_axis)

    assert np.array_equal(mid, expect)


def test_month_bounds():
    time_axis = xr.DataArray(
        xr.cftime_range("2005-01", periods=12, freq="MS"), dims="time"
    )

    expect = np.array(
        [
            [
                cftime.DatetimeGregorian(2005, 1, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 2, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 2, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 3, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 3, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 4, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 4, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 5, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 5, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 6, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 6, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 7, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 7, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 8, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 8, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 9, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 9, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 10, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 10, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 11, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 11, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2005, 12, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
            [
                cftime.DatetimeGregorian(2005, 12, 1, 0, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2006, 1, 1, 0, 0, 0, 0, has_year_zero=False),
            ],
        ],
        dtype=object,
    )

    mid = utils.month_bounds(time_axis)

    assert np.array_equal(mid, expect)


def test_cmorizer_fx():
    ds = cx.cordex_domain("EUR-11", dummy="topo")
    filename = cmor.cmorize_variable(
        ds,
        "orog",
        mapping_table={"orog": {"varname": "topo"}},
        cmor_table=cx.tables.cordex_cmor_table("CORDEX_fx"),
        dataset_table=cx.tables.cordex_cmor_table("CORDEX_remo_example"),
        grids_table=cx.tables.cordex_cmor_table("CORDEX_grids"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert "orog" in output


def test_cmorizer_mon():
    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_mon")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = cmor.cmorize_variable(
        ds,
        "tas",
        mapping_table={"tas": {"varname": "TEMP2"}},
        cmor_table=cx.tables.cordex_cmor_table("CORDEX_mon"),
        dataset_table=cx.tables.cordex_cmor_table("CORDEX_remo_example"),
        grids_table=cx.tables.cordex_cmor_table("CORDEX_grids"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert output.dims["time"] == 12
    assert "tas" in output


@pytest.mark.parametrize("table, tdim", [("CORDEX_day", 3), ("CORDEX_1hr", 49)])
def test_cmorizer_subdaily(table, tdim):
    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_1hr")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = cmor.cmorize_variable(
        ds,
        "tas",
        mapping_table={"tas": {"varname": "TEMP2"}},
        cmor_table=cx.tables.cordex_cmor_table(table),
        dataset_table=cx.tables.cordex_cmor_table("CORDEX_remo_example"),
        grids_table=cx.tables.cordex_cmor_table("CORDEX_grids"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
        allow_resample=True,
    )
    output = xr.open_dataset(filename)
    assert "tas" in output
    assert output.dims["time"] == tdim
