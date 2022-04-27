from pooch import retrieve

cache_url = "~/.py-cordex"


def fetch_vg2500():
    """Fetch Germany Verwaltungsgebiete 1:2,500,000"""
    fname = retrieve(
        path=cache_url,
        url="https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/2020/vg2500_01-01.gk3.shape.zip",
        known_hash="md5:5a1a86cd131decd9cf116dbfc1a66f17",
    )
    return fname


def fetch_prudence():
    """Fetch Prudence regions table"""
    fname = retrieve(
        path=cache_url,
        url="https://raw.githubusercontent.com/euro-cordex/tables/master/regions/prudence.csv",
        known_hash="d87691a873110c9e3e4460a0ed35cd15f11f2a42aa86aced76feae9e87e8bed2",
    )
    return fname
