import datetime as dt
import os

import cftime
import cftime as cfdt
import numpy as np
import pandas as pd
import pytest
import xarray as xr

import cordex as cx
from cordex import cmor
from cordex.cmor import utils
from cordex.tables import cordex_cmor_table
from cordex.cmor.cmor import _crop_to_cordex_domain

from . import requires_pint_xarray

mapping_table = {"orog": {"varname": "topo"}, "tas": {"varname": "TEMP2"}}

table_prefix = "CORDEX-CMIP6"

cmor.set_options(table_prefix=table_prefix)


def create_sdepth_ds():
    ds = cx.domain("EUR-11", dummy="topo")
    ds = ds.drop("topo").assign(
        tsl=ds.topo.expand_dims(
            time=pd.date_range("2000-01-01T12:00:00", periods=3, freq="D"),
            sdepth=[0.05, 0.10, 0.2],
        )
    )
    ds.tsl.attrs["units"] = "K"
    ds.sdepth.attrs["axis"] = "Z"
    ds.sdepth.attrs["units"] = "m"
    return ds


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


@pytest.mark.parametrize("domain_id", ["EUR-11", "SAM-44", "AFR-22"])
def test_crop_to_domain(domain_id):
    ds = cx.domain(domain_id)

    cropped = _crop_to_cordex_domain(ds, domain_id)
    assert ds.equals(cropped)
    # assert larger domain is correctly cropped
    pad = ds.pad(rlon=(1, 1), rlat=(1, 1), mode="reflect", reflect_type="odd")
    cropped = _crop_to_cordex_domain(pad, domain_id)
    assert ds.equals(cropped)

    # check for tolerance
    cropped = _crop_to_cordex_domain(
        pad.reindex(rlon=pad.rlon * 1.00001, method="nearest"), domain_id
    ).assign_coords(rlon=ds.rlon, rlat=ds.rlat)
    assert ds.equals(cropped)


def run_cmorizer(ds, out_name, domain_id, table_id, dataset_table=None, **kwargs):
    if dataset_table is None:
        dataset_table = cordex_cmor_table(f"{table_prefix}_remo_example")
    return cmor.cmorize_variable(
        ds,
        out_name,
        mapping_table=mapping_table,
        cmor_table=cordex_cmor_table(f"{table_prefix}_{table_id}"),
        dataset_table=dataset_table,
        grids_table=cordex_cmor_table(f"{table_prefix}_grids"),
        domain_id=domain_id,
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
    assert output.dims == {"rlat": 412, "rlon": 424, "vertices": 4}


def test_cmorizer_mon():
    ds = cx.tutorial.open_dataset("remo_EUR-11_TEMP2_mon")
    filename = run_cmorizer(ds, "tas", "EUR-11", "mon")
    output = xr.open_dataset(filename)
    # assert output.dims["time"] == 12
    assert "tas" in output
    assert output.dims == {
        "time": 12,
        "bnds": 2,
        "rlat": 412,
        "rlon": 424,
        "vertices": 4,
    }


def test_cmorizer_mon_sdepth():
    ds = create_sdepth_ds()
    filename = run_cmorizer(ds, "tsl", "EUR-11", "day")
    return filename


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
    table = f"{table_prefix}_{table_id}"
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


def test_cmorizer_with_xwrf():
    import xwrf

    ds_old = xwrf.tutorial.open_dataset("wrfout")
    ds_new = ds_old.xwrf.postprocess()
    # we use the lowest layer as tas dummy
    ds = (
        ds_new[
            ["x", "y", "air_potential_temperature", "XLAT", "XLONG", "wrf_projection"]
        ]
        .isel(z=0)
        .rename(Time="time", air_potential_temperature="tas")
    )
    # extend time axis and make dummy monthly dataset
    time = pd.date_range("2025-01-01", periods=12, freq="MS")
    ds = (
        xr.decode_cf(ds, decode_coords="all")
        .squeeze(drop=True)
        .expand_dims(time=time, axis=0)
    )
    table_id = "mon"
    filename = cmor.cmorize_variable(
        ds,
        "tas",
        cmor_table=cordex_cmor_table(f"{table_prefix}_{table_id}"),
        dataset_table=cordex_cmor_table(f"{table_prefix}_remo_example"),
        grids_table=cordex_cmor_table(f"{table_prefix}_grids"),
        domain_id="West-Coast",
        crop=False,
    )
    ds_out = xr.open_dataset(filename)
    assert "tas" in ds_out
    assert ds_out.dims["time"] == 12

    grid_mapping = ds_out.cf["grid_mapping"]
    assert grid_mapping.grid_mapping_name == "lambert_conformal_conic"
    expected_grid_mapping_attrs = {
        "grid_mapping_name": "lambert_conformal_conic",
        "standard_parallel": np.array([30.0, 60.0]),
        "longitude_of_central_meridian": np.float64(-70.0),
        "latitude_of_projection_origin": np.float64(38.0),
        "false_easting": np.float64(0.0),
        "false_northing": np.float64(0.0),
    }
    for key, value in expected_grid_mapping_attrs.items():
        assert key in grid_mapping.attrs
        if isinstance(value, np.ndarray):
            assert np.array_equal(value, grid_mapping.attrs[key])
        else:
            assert value == grid_mapping.attrs[key]
