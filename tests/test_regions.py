import cordex as cx


# naive download test
def test_germany():
    cx.regions.germany.geodataframe()
    cx.regions.germany.regionmask()


# naive download test
def test_prudence():
    cx.regions.prudence.geodataframe
    cx.regions.prudence.regionmask
