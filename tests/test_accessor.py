import pytest
import xarray as xr
import numpy as np

import cordex as cx
from cordex.accessor import CordexDataArrayAccessor, CordexDatasetAccessor  # noqa


@pytest.mark.parametrize(
    "domain_id", ["EUR-11", "EUR-22", "EUR-44", "EUR-11i", "AFR-44"]
)
def test_guess_info(domain_id):
    ds = xr.decode_cf(cx.cordex_domain(domain_id, dummy=True), decode_coords="all")
    info = cx.domain_info(domain_id)
    assert set(ds.cx.info().values()) <= set(info.values())
    assert set(ds.dummy.cx.info().values()) <= set(info.values())
    assert ds.cx.domain_id == domain_id
    assert ds.cx.guess() == info
    assert ds.dummy.cx.guess() == info

    # assert we guess unknown domains
    ds = xr.decode_cf(cx.cordex_domain(domain_id, dummy=True), decode_coords="all")
    del ds.attrs["CORDEX_domain"]
    assert ds.cx.domain_id == domain_id


@pytest.mark.parametrize(
    "filename",
    [
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1_mon_197901-198012",
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012",
    ],
)
def test_dataset_guess(filename):
    ds = cx.tutorial.open_dataset(filename)
    assert ds.cx.guess() == cx.domain_info(ds.cx.domain_id)


@pytest.mark.parametrize("domain_id", ["EUR-11", "EUR-44", "SAM-44", "AFR-22"])
def test_rewrite_coords(domain_id):
    # Create a sample dataset
    grid = cx.domain(domain_id)

    # Create typical coordinate precision issue by adding random noise
    rlon_noise = np.random.randn(*grid.rlon.shape) * np.finfo("float32").eps
    rlat_noise = np.random.randn(*grid.rlat.shape) * np.finfo("float32").eps
    # lon_noise = np.random.randn(*grid.lon.shape) * np.finfo("float32").eps
    # lat_noise = np.random.randn(*grid.lat.shape) * np.finfo("float32").eps

    # Call the rewrite_coords function for "xy" coordinates
    grid_noise = grid.assign_coords(
        rlon=grid.rlon + rlon_noise, rlat=grid.rlat + rlat_noise
    )
    rewritten_data = grid_noise.cx.rewrite_coords(coords="xy")

    np.testing.assert_array_equal(rewritten_data.rlon, grid.rlon)
    np.testing.assert_array_equal(rewritten_data.rlat, grid.rlat)
    xr.testing.assert_identical(rewritten_data, grid)
