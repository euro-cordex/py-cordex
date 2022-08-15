import datetime as dt

import cftime as cfdt
import pytest
import xarray as xr

import cordex as cx
from cordex import cmor


def test_cftime():
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


def test_cmorizer_fx():
    ds = cx.cordex_domain("EUR-11", dummy="topo")
    filename = cmor.cmorize_variable(
        ds,
        "orog",
        mapping_table={"orog": {"varname": "topo"}},
        cmor_table=cx.tables.cmip6_cmor_table("CMIP6_fx"),
        dataset_table=cx.tables.cordex_cmor_table("CORDEX_remo_example"),
        grids_table=cx.tables.cmip6_cmor_table("CMIP6_grids"),
        CORDEX_domain="EUR-11",
        time_units=None,
        allow_units_convert=True,
    )
    output = xr.open_dataset(filename)
    assert "orog" in output
