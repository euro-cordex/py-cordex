import datetime as dt
import os

import cftime
import cftime as cfdt
import numpy as np
import pytest
import xarray as xr

import cordex as cx
from cordex import cmor
from cordex.cmor import utils
from cordex.tables import cordex_cmor_table

from . import requires_pint_xarray

mapping_table = {"orog": {"varname": "topo"}, "tas": {"varname": "TEMP2"}}


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


def run_cmorizer(ds, out_name, domain_id, table_id, dataset_table=None, **kwargs):
    if dataset_table is None:
        dataset_table = cordex_cmor_table("CORDEX_remo_example")
    return cmor.cmorize_variable(
        ds,
        out_name,
        mapping_table=mapping_table,
        cmor_table=cordex_cmor_table(f"CORDEX_{table_id}"),
        dataset_table=dataset_table,
        grids_table=cordex_cmor_table("CORDEX_grids"),
        CORDEX_domain=domain_id,
        replace_coords=True,
        allow_units_convert=True,
        allow_resample=True,
        **kwargs,
    )


def test_cmorizer_fx():
    ds = cx.cordex_domain("EUR-11", dummy="topo")
    filename = run_cmorizer(ds, "orog", "EUR-11", "fx")
    output = xr.open_dataset(filename)
    assert "orog" in output


def test_cmorizer_mon():
    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_mon")
    filename = run_cmorizer(ds, "tas", "EUR-11", "mon")
    output = xr.open_dataset(filename)
    assert output.dims["time"] == 12
    assert "tas" in output


@pytest.mark.parametrize("table_id, tdim", [("day", 3), ("1hr", 49)])
def test_cmorizer_subdaily(table_id, tdim):
    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_1hr")
    eur11 = cx.cordex_domain("EUR-11")
    ds = ds.assign_coords({"lon": eur11.lon, "lat": eur11.lat})
    filename = run_cmorizer(ds, "tas", "EUR-11", table_id)
    output = xr.open_dataset(filename)
    assert "tas" in output
    assert output.dims["time"] == tdim


@requires_pint_xarray
def test_units_convert():
    import pint_xarray  # noqa
    from cf_xarray.units import units  # noqa

    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_mon")
    ds["TEMP2"] = ds.TEMP2.pint.quantify().pint.to("degC").pint.dequantify(format="cf")
    filename = run_cmorizer(ds, "tas", "EUR-11", "mon")
    output = xr.open_dataset(filename)
    assert output.tas.units == "K"


@pytest.mark.parametrize("table_id", ["1hr", "6hr", "day", "mon"])
def test_table_id(table_id):
    table = f"CORDEX_{table_id}"
    filename = cordex_cmor_table(table)
    tid = utils.get_table_id(utils._read_table(filename))
    assert tid == table_id


def test_table_manipulation():
    home = os.environ["HOME"]
    print(f"writing to {home}")
    ds = cx.cordex_domain("EUR-11", dummy="topo")
    filename = run_cmorizer(ds, "orog", "EUR-11", "fx", outpath=home)
    print(filename)
    assert home in filename
