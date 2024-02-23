import numpy as np
import pytest
import xarray as xr

import cordex as cx

from . import requires_cartopy

# from cordex.utils import _get_info, _guess_domain


@pytest.mark.parametrize("domain_id", ["EUR-11", "EUR-11i", "SAM-44", "AFR-22"])
@pytest.mark.parametrize("bounds", [False, True])
@pytest.mark.parametrize("mapping_name", [None, "rotated_pole"])
@pytest.mark.parametrize("dummy", [False, True, "data", "topo"])
@pytest.mark.parametrize("cell_area", [False, True])
def test_domain_coordinates(domain_id, bounds, mapping_name, dummy, cell_area):
    ds = cx.cordex_domain(
        domain_id,
        bounds=bounds,
        mapping_name=mapping_name,
        dummy=dummy,
        cell_area=cell_area,
    )
    assert ds.cf["X"].ndim == 1
    assert ds.cf["Y"].ndim == 1

    if ds.cf.grid_mapping_names:
        assert ds.cf["longitude"].ndim == 2
        assert ds.cf["latitude"].ndim == 2
        if mapping_name:
            assert [mapping_name] in ds.cf.grid_mapping_names.values()
        else:
            assert ["rotated_latitude_longitude"] in ds.cf.grid_mapping_names.values()
    else:
        assert ds.cf["longitude"].ndim == 1
        assert ds.cf["latitude"].ndim == 1
        assert ds.cf.grid_mapping_names == {}

    if bounds is True:
        assert "longitude" in ds.cf.bounds
        assert "latitude" in ds.cf.bounds

    if cell_area is True:
        assert "areacella" in ds

    if dummy is True:
        assert "dummy" in ds
    if isinstance(dummy, str):
        assert dummy in ds


@pytest.mark.parametrize("domain_id", ["EUR-11", "EUR-11i", "SAM-44", "AFR-22"])
def test_domain(domain_id):
    ds = cx.cordex_domain(domain_id)
    ds.attrs["CORDEX_domain"] == domain_id
    # test attributes
    assert "institution" in cx.cordex_domain(domain_id, attrs="CORDEX").attrs
    # ensure rounding errors fixed


def test_constructor():
    eur11 = cx.cordex_domain("EUR-11")
    eur11_user = cx.create_dataset(
        nlon=424,
        nlat=412,
        dlon=0.11,
        dlat=0.11,
        ll_lon=-28.375,
        ll_lat=-23.375,
        pollon=-162.00,
        pollat=39.25,
    )
    assert eur11_user.equals(eur11)
    assert np.float64(eur11.rlon.isel(rlon=34).data) == -24.635


@pytest.mark.parametrize("bounds", [False, True])
# @pytest.mark.parametrize("", [2, 3])
def test_domain_info(bounds):
    import pandas as pd

    info = {
        "short_name": "EUR-11",
        "region": 4,
        "long_name": "Europe",
        "nlon": 424,
        "nlat": 412,
        "ll_lon": -28.375,
        "ur_lon": 18.155,
        "ll_lat": -23.375,
        "ur_lat": 21.835,
        "dlon": 0.11,
        "dlat": 0.11,
        "pollon": -162.0,
        "pollat": 39.25,
    }
    table = pd.DataFrame(info, index=[0]).set_index("short_name")
    assert cx.domain_info("EUR-11", table) == info
    xr.testing.assert_equal(
        cx.create_dataset(**info, bounds=bounds),
        cx.cordex_domain("EUR-11", bounds=bounds),
    )


@requires_cartopy
def test_vertices():
    eur11 = cx.cordex_domain("EUR-11")
    import cartopy.crs as ccrs

    pole = (
        eur11.rotated_latitude_longitude.grid_north_pole_longitude,
        eur11.rotated_latitude_longitude.grid_north_pole_latitude,
    )
    cx.vertices(eur11.rlon, eur11.rlat, src_crs=ccrs.RotatedPole(*pole))
