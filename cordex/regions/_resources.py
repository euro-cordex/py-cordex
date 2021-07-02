from pooch import retrieve
import pandas
import pooch

cache_url = "~/.py-cordex"

# Create a new friend to manage your sample data storage
REGION_RESOURCE = pooch.create(
    # Folder where the data will be stored. For a sensible default, use the
    # default cache folder for your OS.
    path=cache_url, #pooch.os_cache("py-cordex"),
    # Base URL of the remote data store. Will call .format on this string
    # to insert the version (see below).
    base_url="https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/",
    # Pooches are versioned so that you can use multiple versions of a
    # package simultaneously. Use PEP440 compliant version number. The
    # version will be appended to the path.
    # version=version,
    # If a version as a "+XX.XXXXX" suffix, we'll assume that this is a dev
    # version and replace the version with this string.
    version_dev="master",
    # An environment variable that overwrites the path.
    env="PY_CORDEX_DATA_DIR",
    # The cache file registry. A dictionary with all files managed by this
    # pooch. Keys are the file names (relative to *base_url*) and values
    # are their respective SHA256 hashes. Files will be downloaded
    # automatically when needed (see fetch_gravity_data).
    registry={"vg2500_01-01.gk3.shape.zip": "md5:5a1a86cd131decd9cf116dbfc1a66f17"},
)
# You can also load the registry from a file. Each line contains a file
# name and it's sha256 hash separated by a space. This makes it easier to
# manage large numbers of data files. The registry file should be packaged
# and distributed with your software.
# GOODBOY.load_registry(
#    pkg_resources.resource_stream("mypackage", "registry.txt")
# )


# Define functions that your users can call to get back the data in memory
def fetch(dataset):
    """
    Load dataset.
    """
    # Fetch the path to a file in the local storage. If it's not there,
    # we'll download it.
    fname = GOODBOY.fetch(dataset)
    # Load it with numpy/pandas/etc
    # data = pandas.read_csv(fname)
    return fname


def fetch_vg2500():
    """Fetch Germany Verwaltungsgebiete 1:2,500,000"""
    fname = retrieve(
        path=cache_url,
        url="https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/vg2500_01-01.gk3.shape.zip",
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
