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
