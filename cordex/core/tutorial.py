"""
Useful for:

* users learning py-cordex

"""
### code stolen from xarray, I am sorry!
import os
import pathlib

import numpy as np

from xarray import open_dataset as _open_dataset
from xarray import DataArray
from xarray import Dataset

_default_cache_dir_name = "py-cordex_tutorial_data"
base_url = "https://github.com/euro-cordex/py-cordex-data"
version = "main"


external_urls = {}  # type: dict
external_rasterio_urls = {
    "RGB.byte": "https://github.com/mapbox/rasterio/raw/1.2.1/tests/data/RGB.byte.tif",
    "shade": "https://github.com/mapbox/rasterio/raw/1.2.1/tests/data/shade.tif",
}
file_formats = {
    "air_temperature": 3,
    "rasm": 3,
    "ROMS_example": 4,
    "tiny": 3,
    "eraint_uvz": 3,
}


def _construct_cache_dir(path):
    import pooch

    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)

    return path



# idea borrowed from Seaborn
def open_dataset(
    name,
    cache=True,
    cache_dir=None,
    *,
    engine=None,
    **kws,
):
    """
    Open a dataset from the online repository (requires internet).
    If a local copy is found then always use that to avoid network traffic.
    Available datasets:
    * ``"air_temperature"``: NCEP reanalysis subset
    * ``"rasm"``: Output of the Regional Arctic System Model (RASM)
    * ``"ROMS_example"``: Regional Ocean Model System (ROMS) output
    * ``"tiny"``: small synthetic dataset with a 1D data variable
    * ``"era5-2mt-2019-03-uk.grib"``: ERA5 temperature data over the UK
    * ``"eraint_uvz"``: data from ERA-Interim reanalysis, monthly averages of upper level data
    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'air_temperature'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : bool, optional
        If True, then cache data locally for use on subsequent calls
    **kws : dict, optional
        Passed to xarray.open_dataset
    See Also
    --------
    xarray.open_dataset
    """
    try:
        import pooch
    except ImportError as e:
        raise ImportError(
            "tutorial.open_dataset depends on pooch to download and manage datasets."
            " To proceed please install pooch."
        ) from e

    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    cache_dir = _construct_cache_dir(cache_dir)
    if name in external_urls:
        url = external_urls[name]
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"

        url = f"{base_url}/raw/{version}/{path.name}"

    # retrieve the file
    filepath = pooch.retrieve(url=url, known_hash=None, path=cache_dir)
    ds = _open_dataset(filepath, engine=engine, **kws)
    if not cache:
        ds = ds.load()
        pathlib.Path(filepath).unlink()

    return 


# idea borrowed from Seaborn
def open_dataset(
    name,
    cache=True,
    cache_dir=None,
    *,
    engine=None,
    **kws,
):
    """
    Open a dataset from the online repository (requires internet).
    If a local copy is found then always use that to avoid network traffic.
    Available datasets:
    * ``"air_temperature"``: NCEP reanalysis subset
    * ``"rasm"``: Output of the Regional Arctic System Model (RASM)
    * ``"ROMS_example"``: Regional Ocean Model System (ROMS) output
    * ``"tiny"``: small synthetic dataset with a 1D data variable
    * ``"era5-2mt-2019-03-uk.grib"``: ERA5 temperature data over the UK
    * ``"eraint_uvz"``: data from ERA-Interim reanalysis, monthly averages of upper level data
    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'air_temperature'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : bool, optional
        If True, then cache data locally for use on subsequent calls
    **kws : dict, optional
        Passed to xarray.open_dataset
    See Also
    --------
    xarray.open_dataset
    """
    try:
        import pooch
    except ImportError as e:
        raise ImportError(
            "tutorial.open_dataset depends on pooch to download and manage datasets."
            " To proceed please install pooch."
        ) from e

    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    cache_dir = _construct_cache_dir(cache_dir)
    if name in external_urls:
        url = external_urls[name]
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"

        url = f"{base_url}/raw/{version}/{path.name}"

    # retrieve the file
    filepath = pooch.retrieve(url=url, known_hash=None, path=cache_dir)
    ds = _open_dataset(filepath, engine=engine, **kws)
    if not cache:
        ds = ds.load()
        pathlib.Path(filepath).unlink()

    return ds
