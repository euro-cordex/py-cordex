import pandas as pd
import pooch


base_url = "https://raw.githubusercontent.com/euro-cordex/tables/master/"

cache_url = "~/.py-cordex"

DOMAIN_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "domains/",
    registry={
        "cordex.csv": "115acce1612ff7d0908ba6693d3ef3892ee01f8656239b7383abdd633a3cbd7d",
        "cordex-high-res.csv": "a9b69a4a0db7e0b5920fc5ae4b030b770697f33e993f4b8fe39bd3a8093bad7c",
        "cordex-fps.csv": "09185a260a505c64edf1970a0c7ae91f5adb4d6f08d7e7af83767b96bd5017e3",
        "cordex-core.csv": "fd81e159ae5540581424c2f43a462b324e37b3a2cd7deddc2cf415c077fd8c2b",
    },
)


REGION_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "regions/",
    registry={
        "prudence.csv": "d87691a873110c9e3e4460a0ed35cd15f11f2a42aa86aced76feae9e87e8bed2",
    },
)


ECMWF_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "ecmwf/",
    registry={
        "ecmwf_128.csv": "7c29b13c2d74b9442b532268f7f4ed9eaf507cc21f4dc78242495bb387ce6ff4",
    },
)


def fetch_remote_table(name, resource):
    """
    uses pooch to cache files
    """

    # the file will be downloaded automatically the first time this is run.
    return resource.fetch(name)


def read_remote_table(name, resource, index_col=None):

    fname = fetch_remote_table(name, resource)

    return pd.read_csv(fname, index_col=index_col)


def read_region_table(name):

    return read_remote_table(name, resource=REGION_RESOURCE, index_col="area")


def read_cordex_domain_tables():

    resource = DOMAIN_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(table, resource, index_col="short_name")
        for table in resource.registry.keys()
    }


def region_tables():

    return {
        table.split(".")[0]: read_region_table(table)
        for table in REGION_RESOURCE.registry.keys()
    }


def ecmwf_tables():

    resource = ECMWF_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(table, resource, index_col="code")
        for table in ECMWF_RESOURCE.registry.keys()
    }
