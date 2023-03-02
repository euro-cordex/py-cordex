import pytest
import xarray as xr

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
