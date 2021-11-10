import pytest
import cordex as cx


# naive download test
def test_germany():
    mask = cx.regions.germany.geodataframe()
    mask = cx.regions.germany.regionmask()


# naive download test
def test_prudence():
    mask = cx.regions.prudence.geodataframe
    mask = cx.regions.prudence.regionmask
