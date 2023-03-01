import os
from pathlib import Path

import pandas as pd
import pooch

base_url = "https://raw.githubusercontent.com/euro-cordex/tables/main/"

cache_url = "~/.py-cordex"

_default_cache_dir_name = "py-cordex-tables"


def _construct_cache_dir(path):
    import pooch

    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)

    return path


DOMAIN_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "domains/",
    registry={
        "cordex.csv": None,
        "cordex-high-res.csv": None,
        "cordex-fps.csv": None,
        "cordex-core.csv": None,
        "cordex-regular.csv": None,
    },
)


REGION_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "regions/",
    registry={
        "prudence.csv": None,
    },
)


ECMWF_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "ecmwf/",
    registry={
        "ecmwf_128.csv": None,
    },
)


cmor_tables_inpath = str(pooch.os_cache("cmor-tables"))


def fetch_cordex_cmor_table(table):
    return retrieve_cmor_table(
        # table, url="https://github.com/euro-cordex/cordex-cmor-tables/raw/main/Tables"
        table,
        url="https://github.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/raw/main/Tables",
    )


def fetch_cmip6_cmor_table(table):
    return retrieve_cmor_table(
        table, url="https://github.com/PCMDI/cmip6-cmor-tables/raw/master/Tables"
    )


def retrieve_cmor_table(table, url):
    path = cmor_tables_inpath
    if Path(table).suffix == "":
        fname = "{}.json".format(table)
    else:
        fname = table
    return pooch.retrieve(
        os.path.join(url, fname), known_hash=None, fname=fname, path=path
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
