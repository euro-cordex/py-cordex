import numpy as np
import pytest

import cordex as cx
from cordex.utils import cell_area


@pytest.mark.parametrize("domain_id", ["EUR-11", "EUR-11i", "SAM-44", "AFR-22"])
def test_cell_area(domain_id):
    """compare against cdo"""
    from cdo import Cdo

    cdo = Cdo()

    ds = cx.cordex_domain(domain_id, dummy=True, bounds=True)
    area_cx = cell_area(ds)
    area_cdo = cdo.gridarea(input=ds, returnXDataset=True)["cell_area"]

    np.testing.assert_allclose(area_cx, area_cdo, rtol=1.0e-4, atol=0)
