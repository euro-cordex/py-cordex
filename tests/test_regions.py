import pytest
import requests
import cordex as cx

from . import requires_geopandas


SERVICE_URL = "https://daten.gdz.bkg.bund.de"


def _service_available(url: str = SERVICE_URL, timeout: float = 2.0) -> bool:
    """Return True if remote service can be reached quickly.

    Uses a HEAD request first (cheap); falls back to GET if HEAD not allowed.
    Swallows network exceptions and returns False on failure.
    """
    try:
        try:
            r = requests.head(url, timeout=timeout, allow_redirects=True)
            if r.ok:
                return True
            # some servers may not implement HEAD correctly
        except requests.RequestException:
            pass
        r = requests.get(url, timeout=timeout)
        return r.ok
    except requests.RequestException:
        return False


# naive download test
@requires_geopandas
def test_germany():
    if not _service_available():
        pytest.skip("Service unavailable")
    cx.regions.germany.geodataframe()
    cx.regions.germany.regionmask()


# naive download test
@requires_geopandas
def test_prudence():
    cx.regions.prudence.geodataframe
    cx.regions.prudence.regionmask
