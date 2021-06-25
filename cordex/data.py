from pooch import retrieve

##
##
### Download the file and save it locally. Running this again will not cause
### a download. Pooch will check the hash (checksum) of the downloaded file
### against the given value to make sure it's the right file (not corrupted
### or outdated).
##fname = retrieve(
##    url="https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/vg2500_01-01.gk3.shape.zip",
##    known_hash="md5:5a1a86cd131decd9cf116dbfc1a66f17",
##)


"""
Module mypackage/datasets.py
"""
import pkg_resources
import pandas
import pooch

# Get the version string from your project. You have one of these, right?
from . import __version__ as version


# Create a new friend to manage your sample data storage
GOODBOY = pooch.create(
    # Folder where the data will be stored. For a sensible default, use the
    # default cache folder for your OS.
    path=pooch.os_cache("py-cordex"),
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
        url="https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/vg2500_01-01.gk3.shape.zip",
        known_hash="md5:5a1a86cd131decd9cf116dbfc1a66f17",
    )
    return fname
