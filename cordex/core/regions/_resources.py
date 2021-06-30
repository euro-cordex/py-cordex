import pandas as pd
import pooch


base_url = "https://raw.githubusercontent.com/euro-cordex/tables/master/"

cache_url = "~/.cordex"


REGION_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url, #pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url= base_url + "regions/",
    registry={
        "prudence.csv": "d87691a873110c9e3e4460a0ed35cd15f11f2a42aa86aced76feae9e87e8bed2",       
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


def region_tables():
    
    return {table: read_region_table(table) for table in REGION_RESOURCE.registry.keys()}
        
