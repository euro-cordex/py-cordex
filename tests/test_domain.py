import numpy as np
import pytest
import xarray as xr
import pandas as pd

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


@pytest.mark.parametrize("domain_id", ["EUR-11", "EUR-44", "SAM-44", "AFR-22"])
def test_rewrite_coords(domain_id):
    """
    Test the rewrite_coords function for different domains.

    This function tests the rewrite_coords function by creating a sample dataset
    with typical coordinate precision issues (random noise) and verifying that
    the coordinates are correctly rewritten.

    Parameters
    ----------
    domain_id : str
        The domain identifier used to obtain grid information for testing.
    """
    # Create a sample dataset
    grid = cx.domain(domain_id)

    # Create typical coordinate precision issue by adding random noise
    rlon_noise = np.random.randn(*grid.rlon.shape) * np.finfo("float32").eps
    rlat_noise = np.random.randn(*grid.rlat.shape) * np.finfo("float32").eps
    lon_noise = np.random.randn(*grid.lon.shape) * np.finfo("float32").eps
    lat_noise = np.random.randn(*grid.lat.shape) * np.finfo("float32").eps

    # Call the rewrite_coords function for "xy" coordinates
    grid_noise = grid.assign_coords(
        rlon=grid.rlon + rlon_noise, rlat=grid.rlat + rlat_noise
    )
    rewritten_data = cx.rewrite_coords(grid_noise, coords="xy")

    np.testing.assert_array_equal(rewritten_data.rlon, grid.rlon)
    np.testing.assert_array_equal(rewritten_data.rlat, grid.rlat)
    xr.testing.assert_identical(rewritten_data, grid)

    # Call the rewrite_coords function for "lonlat" coordinates
    grid_noise = grid.assign_coords(lon=grid.lon + lon_noise, lat=grid.lat + lat_noise)
    rewritten_data = cx.rewrite_coords(grid_noise, coords="lonlat")

    np.testing.assert_array_equal(rewritten_data.lon, grid.lon)
    np.testing.assert_array_equal(rewritten_data.lat, grid.lat)
    xr.testing.assert_identical(rewritten_data, grid)

    # Call the rewrite_coords function for "all" coordinates
    grid_noise = grid.assign_coords(
        rlon=grid.rlon + rlon_noise,
        rlat=grid.rlat + rlat_noise,
        lon=grid.lon + lon_noise,
        lat=grid.lat + lat_noise,
    )
    rewritten_data = cx.rewrite_coords(grid_noise, coords="all")

    np.testing.assert_array_equal(rewritten_data.rlon, grid.rlon)
    np.testing.assert_array_equal(rewritten_data.rlat, grid.rlat)
    np.testing.assert_array_equal(rewritten_data.lon, grid.lon)
    np.testing.assert_array_equal(rewritten_data.lat, grid.lat)
    xr.testing.assert_identical(rewritten_data, grid)

    grid = cx.domain(domain_id, bounds=True, mip_era="CMIP6")
    grid["vertices_lon"][:] = 0.0
    grid["vertices_lat"][:] = 0.0
    grid.vertices_lon.attrs["hello"] = "world"
    grid.vertices_lat.attrs["hello"] = "world"

    rewritten_data = cx.rewrite_coords(grid, bounds=True)
    grid = cx.domain(domain_id, bounds=True, mip_era="CMIP6")

    np.testing.assert_array_equal(rewritten_data.vertices_lon, grid.vertices_lon)
    np.testing.assert_array_equal(rewritten_data.vertices_lat, grid.vertices_lat)

    # check if attributes are now overwritten
    assert rewritten_data.vertices_lon.attrs["hello"] == "world"
    assert rewritten_data.vertices_lat.attrs["hello"] == "world"
