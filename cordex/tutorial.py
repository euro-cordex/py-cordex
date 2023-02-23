"""
Useful for:

* users learning py-cordex

"""
# code stolen from xarray, I am sorry!
import os
import pathlib

from xarray import open_dataset as _open_dataset

from .preprocessing import cordex_dataset_id

_default_cache_dir_name = "py-cordex_tutorial_data"
base_url = "https://github.com/euro-cordex/py-cordex-data"
version = "main"


external_urls = {}  # type: dict

file_formats = {
    "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012.nc": 3,
    "rasm": 3,
    "ROMS_example": 4,
    "tiny": 3,
    "eraint_uvz": 3,
}


def _check_netcdf_engine_installed(name):
    version = file_formats.get(name)
    if version == 3:
        try:
            import scipy  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either scipy or "
                    "netCDF4 to be installed."
                )
    if version == 4:
        try:
            import h5netcdf  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either h5netcdf "
                    "or netCDF4 to be installed."
                )


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
    * ``"remo_EUR-11_TEMP2_1hr"``: Remo hourly output
    * ``"remo_EUR-11_TEMP2_mon"``: Remo monthly output
    * ``"remo_EUR-44.nc"``: Remo 3D output on model levels
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012"``: Remo output (rotated pole)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_DHMZ-RegCM4-2_v1_mon_198901-199012"``: RegCM4 output (lambert conformal)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-ALADIN53_v1_mon_197901-198012"``: Aladin Output (lambert conformal)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1_mon_197901-198012"``: Racmo Output (rotated pole)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_RMIB-UGent-ALARO-0_v1_mon_198001-198012"``: Alaro output (lambert conformal)

    Parameters
    ----------
    name : str
        Name of the file containing the dataset.
        e.g. 'tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012'
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

        # from pooch import HTTPDownloader
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


def ensemble():
    """Retrieve a mini CORDEX test ensemble.

    Available datasets:

    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012"``: Remo output (rotated pole)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_DHMZ-RegCM4-2_v1_mon_198901-199012"``: RegCM4 output (lambert conformal)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-ALADIN53_v1_mon_197901-198012"``: Aladin Output (lambert conformal)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1_mon_197901-198012"``: Racmo Output (rotated pole)
    * ``"tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_RMIB-UGent-ALARO-0_v1_mon_198001-198012"``: Alaro output (lambert conformal)

    """
    files = [
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_GERICS-REMO2015_v1_mon_197902-198012",
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_DHMZ-RegCM4-2_v1_mon_198901-199012",
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_CNRM-ALADIN53_v1_mon_197901-198012",
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1_mon_197901-198012",
        "tas_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_RMIB-UGent-ALARO-0_v1_mon_198001-198012",
    ]
    return {cordex_dataset_id(ds): ds for ds in [open_dataset(f) for f in files]}
