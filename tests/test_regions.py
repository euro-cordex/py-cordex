import pytest
import requests
import cordex as cx

from . import requires_geopandas


@pytest.mark.skipif(
    not requests.get("https://daten.gdz.bkg.bund.de").ok, reason="Service unavailable"
)
# naive download test
@requires_geopandas
def test_germany():
    cx.regions.germany.geodataframe()
    cx.regions.germany.regionmask()


# naive download test
@requires_geopandas
def test_prudence():
    cx.regions.prudence.geodataframe
    cx.regions.prudence.regionmask
