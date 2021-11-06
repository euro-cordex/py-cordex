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


CMOR_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path="~/.remo_cmor",
    # The remote data is on Github
    base_url="https://raw.githubusercontent.com/ludwiglierhammer/cmor-tables/main/tables/cordex-cmor-tables-test/Tables/",
    registry={
        "CORDEX_Amon.json": "51a16398ee9c1bcd8f71ee709426a3f06048c7a82f77e5dfba766b40248d4f6a",
        "CORDEX_CV.json": "83c5d19ccf5670c9850d585d89802319f52367aca109ec05221fd5df68544a97",
        "CORDEX_coordinate.json": "bf31a847cdad344b124734a5dbcb28dca740bfe496e2f85ee8af654acd213d8e",
        "CORDEX_formula_terms.json": "6f4e7c60b6089cbc873db9a2e158982b83878780ffbc8d8abe9f172c22756023",
        "CORDEX_grids.json": "970bdb5069598be9b422f9522715d97986c5d406970eb3914c177195224bdc5f",
        "CORDEX_remo_example.json": "f2434504cc9e4c438bbe5b11a1f5a63286a4c356036bad08b23ff255dccbd7d0",
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
