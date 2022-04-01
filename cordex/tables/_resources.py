import pandas as pd
import pooch


base_url = "https://raw.githubusercontent.com/euro-cordex/tables/main/"

cache_url = "~/.py-cordex"

cmor_table_version = "0.1.1"

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


CMOR_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path="~/.cordex-cmor-tables",
    # The remote data is on Github
    base_url="https://raw.githubusercontent.com/ludwiglierhammer/cmor-tables/v{}/tables/cordex-cmor-tables-test/Tables/".format(
        cmor_table_version
    ),
    registry={
        "CORDEX_Amon.json": None,
        "CORDEX_day.json": None,
        "CORDEX_1hr.json": None,
        "CORDEX_3hr.json": None,
        "CORDEX_fx.json": None,
        "CORDEX_CV.json": None,
        "CORDEX_coordinate.json": None,
        "CORDEX_formula_terms.json": None,
        "CORDEX_grids.json": None,
        "CORDEX_remo_example.json": None,
    },
)


def fetch_cordex_cmor_table(table):
    fmt = "CORDEX_{}.json"
    return CMOR_RESOURCE.fetch(fmt.format(table))


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
