"""Preprocessing for Cordex models

based on https://github.com/jbusecke/cmip6_preprocessing/blob/master/cmip6_preprocessing/preprocessing.py

"""


import numpy as np
import xarray as xr

from ..domain import cordex_domain
from . import known_issues as fixes

regridder = None


def _init_regridder(src_grid, trg_grid, method="bilinear", **kwargs):
    import xesmf as xe

    global regridder
    regridder = xe.Regridder(src_grid, trg_grid, method=method, **kwargs)


def _maybe_make_list(item):
    "utility function to make sure output is a list"
    if isinstance(item, str):
        return [item]
    elif isinstance(item, list):
        return item
    else:
        return list(item)


def cordex_renaming_dict():
    """a universal renaming dict. Keys correspond to source id (model name)
    and valuse are a dict of target name (key) and a list of variables that
    should be renamed into the target."""
    rename_dict = {
        # dim labels (order represents the priority when checking for the dim labels)
        "lon": ["longitude"],
        "lat": ["latitude"],
        "lev": ["deptht", "olevel", "zlev", "olev", "depth"],
        "bnds": ["bnds", "axis_nbounds", "d2"],
        "lon_vertices": ["longitude_vertices", "lon_bounds"],
        "lat_vertices": ["latitude_vertices", "lat_bounds"],
        "rotated_latitude_longitude": ["rotated_pole"],
        # coordinate labels
        #   "lon": ["longitude", "nav_lon"],
        #   "lat": ["latitude", "nav_lat"],
        "lev_bounds": [
            "deptht_bounds",
            "lev_bnds",
            "olevel_bounds",
            "zlev_bnds",
        ],
        "lon_bounds": [
            "bounds_lon",
            "bounds_nav_lon",
            "lon_bnds",
            "x_bnds",
            "vertices_longitude",
        ],
        "lat_bounds": [
            "bounds_lat",
            "bounds_nav_lat",
            "lat_bnds",
            "y_bnds",
            "vertices_latitude",
        ],
        "time_bounds": ["time_bnds"],
    }
    return rename_dict


def _invert_dict(rdict):
    exploded_dict = {}
    # there is probably a more effective way to 'invert' a dictionary
    for k, v in rdict.items():
        v = _maybe_make_list(v)
        for vv in v:
            exploded_dict[vv] = k
    return exploded_dict


def rename_cordex(ds, rename_dict=None):
    """Homogenizes cordex datasets to common naming conventions.

    Renames dataset names and attributes according to CORDEX conventions.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    rename_dict: dict
        Renaming dictionary, defaults to None. If rename_dict is None,
        a standard CORDEX dictionary is used.

    Returns
    -------
    ds : xr.Dataset
        Dataset according to CORDEX conventions.

    """
    ds = ds.copy()
    attrs = {k: v for k, v in ds.attrs.items()}

    if rename_dict is None:
        rename_dict = cordex_renaming_dict()

    inverted_rename_dict = _invert_dict(rename_dict)

    ds_reset = ds.reset_coords()

    def _maybe_rename(obj, rdict):
        return obj.rename({kk: vv for kk, vv in rdict.items() if kk in obj.dims})

    # first take care of the dims and reconstruct a clean ds
    ds = xr.Dataset(
        {
            k: _maybe_rename(ds_reset[k], inverted_rename_dict)
            for k in ds_reset.data_vars
        }
    )
    # return {
    #        k: _maybe_rename(ds_reset[k], inverted_rename_dict)
    #        for k in ds_reset.data_vars
    #    }
    # return ds
    # ds_reset = ds_reset.assign_coords()

    # special treatment for 'lon'/'lat' if there is no 'x'/'y' after renaming process
    for di, co in [("rlon", "lon"), ("rlat", "lat")]:
        if di not in ds.dims and co in ds.dims:
            ds = ds.rename({co: di})

    # now rename the variables
    # try and pass here, cause some of the datasets (MIROC) have like 3 times the same info
    # e.g. lev/sigma/zlev...not sure this is the best way to handle this with
    # a silent fail here though...
    for va in ds.data_vars:
        try:
            ds = ds.rename({va: inverted_rename_dict[va]})
        except Exception:
            pass

    # handle WRF where lon lat might be in variable dims
    for va in ds.data_vars:
        try:
            if "lon" in ds[va].dims and "lat" in ds[va].dims:
                ds[va] = ds[va].rename({"lon": "rlon", "lat": "rlat"})
        except Exception:
            pass

    #  re-set lon lat to coordinates
    for coord in ["lat", "lon"]:
        if coord in ds.data_vars:
            ds = ds.set_coords(coord)

    # numpy arrays as netcdf attributes do not work with dict_union
    # in intake_esm...
    for va in ds.data_vars:
        for key, value in ds[va].attrs.items():
            if isinstance(value, np.ndarray):
                ds[va].attrs[key] = list(value)
            if key == "grid_mapping" and value == "rotated_pole":
                ds[va].attrs[key] = "rotated_latitude_longitude"

    # restore attributes
    ds.attrs = attrs

    # numpy arrays as netcdf attributes do not work with dict_union
    # in intake_esm...
    for key, value in ds.attrs.items():
        if isinstance(value, np.ndarray):
            ds.attrs[key] = list(value)

    return ds


def fix_known_issues(ds_dict):
    """Fix known issues in CORDEX datasets"""
    for ds_id, ds in ds_dict.items():
        if ds.model_id == "ALADIN53":
            ds_dict[ds_id] = fixes.ALADIN53(ds)
    return ds_dict


def promote_empty_dims(ds):
    """Convert empty dimensions to actual coordinates.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.

    Returns
    -------
    ds : xr.Dataset
        Dataset with promoted dims.

    """
    ds = ds.copy()
    for di in ds.dims:
        if di not in ds.coords:
            ds = ds.assign_coords({di: ds[di]})
    return ds


def attr_to_coord(ds, attr, expand=True):
    """Promote a dataset attribute to a coordinate.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.

    Returns
    -------
    ds : xr.Dataset
        Dataset with new coordinate.

    """
    ds = ds.copy()
    value = ds.attrs[attr]
    ds = ds.assign_coords({attr: value})
    if expand is True:
        return ds.expand_dims(dim=attr)
    return ds


def check_domain(ds, domain=None):
    """Check if coordinates are conform with archive specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    """
    if domain is None:
        domain = ds.attrs.get("CORDEX_domain", None)
    dm = cordex_domain(domain)
    if "rotated_latitude_longitude" in ds:
        assert ds.rlon.size == dm.rlon.size
        assert ds.rlat.size == dm.rlat.size
        assert np.all(ds.rlon == dm.rlon)
        assert np.all(ds.rlat == dm.rlat)
    return True


def replace_rlon_rlat(ds, domain=None):
    """Replace rlon rlat coordinates with archive specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    Returns
    -------
    ds : xr.Dataset
        Dataset with updated rlon, rlat.

    """
    ds = ds.copy()
    if domain is None:
        domain = ds.attrs.get("CORDEX_domain", None)
    dm = cordex_domain(domain)
    for coord in ["rlon", "rlat"]:
        if coord in ds.coords:
            ds = ds.drop(coord)  #
        ds = ds.assign_coords({coord: dm[coord]})
    return ds


def replace_vertices(ds, domain=None):
    """Replace vertices coordinates with archive specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    Returns
    -------
    ds : xr.Dataset
        Dataset with updated vertices.

    """
    ds = ds.copy()
    if domain is None:
        domain = ds.attrs.get("CORDEX_domain", None)
    dm = cordex_domain(domain, add_vertices=True)
    for var in ["lon_vertices", "lat_vertices"]:
        if var in ds.coords:
            ds = ds.drop(var)  #
        ds[var] = dm[var]
    return ds


def replace_lon_lat(ds, domain=None):
    """Replace lon lat coordinates with archive specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    Returns
    -------
    ds : xr.Dataset
        Dataset with updated lon, lat.

    """
    ds = ds.copy()
    if domain is None:
        domain = ds.attrs.get("CORDEX_domain", None)
    dm = cordex_domain(domain)
    for coord in ["lon", "lat"]:
        if coord in ds.coords:
            ds = ds.drop(coord)  #
        ds = ds.assign_coords({coord: dm[coord]})
    return ds


def replace_coords(ds, domain=None):
    """Replace spatial coordinates with coordinates according
    to archive specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    Returns
    -------
    ds : xr.Dataset
        Dataset with updated rlon, rlat, lon and lat.


    """
    ds = ds.copy()
    ds = replace_rlon_rlat(ds, domain)
    ds = replace_lon_lat(ds, domain)
    return ds


def replace_grid(ds, domain=None):
    """Replace complete grid definitions with CORDEX archive
    specifications.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    domain: str
        CORDEX domain name. Required if domain is not specified in the
        dataset attributes (e.g., CORDEX_domain).

    Returns
    -------
    ds : xr.Dataset
        Dataset with updated: rlon, rlat, lon, lat, lon_vertices, lat_vertices.


    """
    ds = ds.copy()
    ds = replace_rlon_rlat(ds, domain)
    ds = replace_lon_lat(ds, domain)
    ds = replace_vertices(ds, domain)
    return ds


def get_grid_mapping_name(ds):
    """Returns grid mapping name.

    Returns the grid_mapping_name attribute of the grid_mapping
    in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.

    Returns
    -------
    grid_mapping_name : str
        grid_mapping_name attribute.

    """
    return get_grid_mapping(ds).grid_mapping_name


def get_grid_mapping_varname(ds):
    """Returns grid mapping dataarray name.

    Returns the variable name of the first grid mapping found
    in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.

    Returns
    -------
    grid_mapping_varname : str
        Variable name of the grid mapping.

    """
    return next(
        ds[va].attrs["grid_mapping"]
        for va in ds.data_vars
        if "grid_mapping" in ds[va].attrs
    )


def get_grid_mapping(ds):
    """Returns grid mapping dataarray.

    Returns the grid_mapping dataarray of the first grid mapping
    found in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.

    Returns
    -------
    grid_mapping : xr.DataArray
        Dataarray containing the grid mapping meta data.

    """
    return ds[get_grid_mapping_varname(ds)]


def remap_lambert_conformal(ds, regridder=None, domain=None):
    """Remap lambert conformal grid to rotated pole grid

    Remap CORDEX dataset with lambert conformal grid mapping to
    rotated pole grid mapping using xesmf.

    Parameters
    ----------
    ds : xr.Dataset
        CORDEX like dataset.
    regridder: xesmf.Regridder
        Regridder instance for use in the remapping. If regridder is None,
        an instance will be created to remap the dataset to the rotated_pole
        grid defined by the CORDEX_domain dataset attribute.
    domain : str
        Domain name of the target grid. If domain is None, remap the dataset
        to the rotated_pole grid defined by the CORDEX_domain dataset attribute.

    Returns
    -------
    remapped dataset : xr.Dataset
        Dataset with lambert conformal data remapped to rotated pole.


    """

    def grid_mapping_name(da):
        try:
            return ds[da.grid_mapping].grid_mapping_name
        except Exception:
            return None

    ds = ds.copy()
    ds_attrs = ds.attrs
    if domain is None:
        domain = ds.attrs.get("CORDEX_domain", None)
    dm = cordex_domain(domain)
    if regridder is None:
        import xesmf as xe

        regridder = xe.Regridder(ds, dm, method="bilinear")
    for va in ds.data_vars:
        try:
            if grid_mapping_name(ds[va]) == "lambert_conformal_conic":
                old_mapping = ds[va].attrs["grid_mapping"]
                attrs = ds[va].attrs
                ds[va] = regridder(ds[va])
                ds[va].attrs.update(attrs)
                ds[va].attrs.update({"grid_mapping": "rotated_latitude_longitude"})
                try:
                    ds = ds.drop(old_mapping)
                    ds["rotated_latitude_longitude"] = dm.rotated_latitude_longitude
                except Exception:
                    pass
        except Exception:
            pass
    ds = replace_coords(ds, domain)
    ds.attrs = ds_attrs
    try:
        ds = ds.drop(("x", "y"))
    except Exception:
        pass
    try:
        ds.attrs["CORDEX_domain"] = dm.attrs["CORDEX_domain"]
    except Exception:
        pass
    return ds


def flatten_coordinate_to_dset_id(ds, coord):
    """Flattens an arbitrary coordinate to become part of the dataset id."""
    flatten = {}
    for xcoord, ds_coord in ds.groupby(coord):
        dset_id = cordex_dataset_id(ds_coord) + ".{}".format(xcoord)
        flatten[dset_id] = ds_coord
    return flatten


def member_id_to_dset_id(ds_dict):
    """Expand the member coordinate into the dataset id

    If there are more than two members in the dataset, the function
    will give back a dict with new dataset ids containing the member
    in the key and a dataset for each member. This is useful for creating
    and ensemble dataset to avoid sparse datasets.

    Parameters
    ----------
    ds_dict : dict
        CORDEX ensemble with dataset ids as keys and datasets as items.

    Returns
    -------
    ds_dict : dict
        CORDEX ensemble with keys expanded to the member id.

    """
    ds_split = {}
    for ds in ds_dict.values():
        ds_split.update(flatten_coordinate_to_dset_id(ds, "member"))
    return ds_split


# def dset_ids_to_coord(ds_dict):
#     """Creates a DataArray from dataset ids"""
#     dset_ids = list(ds_dict.keys())
#     dim = xr.DataArray(
#         dset_ids, dims="dset_id", name="dset_id", coords={"dset_id": dset_ids}
#     )
#     return dim


# def align_time_axis(ds_dict):
#     from datetime import datetime as dt

#     for ds in ds_dict.values():
#         # ds = ds.copy()
#         ds.coords["time"] = [dt(date.year, date.month, 15) for date in ds.time.values]
#     return ds_dict


# def concat_along_dset_id(ds_dict, coords="minimal", compat="override", **kwargs):
#     dset_coord = dset_ids_to_coord(ds_dict)
#     ds_dict = align_time_axis(ds_dict)
#     ds_list = []
#     for ds in ds_dict.values():
#         ds = replace_rlon_rlat(ds)
#         ds = replace_lon_lat(ds)
#         ds_list.append(ds)
#     return xr.concat(ds_list, dim=dset_coord, coords=coords, compat=compat, **kwargs)


def sort_ds_dict_by_attr(ds_dict, attr):
    """Sorts the dataset dict by a certain attribute.

    Creates a dictionary containing the datasets sorted according to
    a dataset attribute (e.g., sort by ``experiment_id``).

    Parameters
    ----------
    ds_dict : dict
        CORDEX ensemble with dataset ids as keys and datasets as items.
    attr : str
        Dataset attribute used for sorting, e.g., ``experiment_id``.

    Returns
    -------
    ds_dict : dict
        CORDEX ensemble sorted by dataset attributes.


    """
    from collections import defaultdict

    dsets_sorted = defaultdict(dict)
    for dset_id, ds in ds_dict.items():
        value = ds.attrs[attr]
        dsets_sorted[value][dset_id] = ds
    return dsets_sorted


def correct_lon(ds):
    """Wraps negative x and lon values around to have 0-360 lons.
    longitude names expected to be corrected with `rename_cordex`"""
    ds = ds.copy()

    # remove out of bounds values found in some
    # models as missing values
    ds["lon"] = ds["lon"].where(abs(ds["lon"]) <= 1000)
    ds["lat"] = ds["lat"].where(abs(ds["lat"]) <= 1000)

    # adjust lon convention
    lon = ds["lon"].where(ds["lon"] > 0, 360 + ds["lon"])
    ds = ds.assign_coords(lon=lon)

    if "lon_bounds" in ds.variables:
        lon_b = ds["lon_bounds"].where(ds["lon_bounds"] > 0, 360 + ds["lon_bounds"])
        ds = ds.assign_coords(lon_bounds=lon_b)

    return ds


def _key_from_attrs(ds, attrs, sep="."):
    return sep.join([ds.attrs[i] if i in ds.attrs.keys() else "none" for i in attrs])


def get_rotated_pole_datasets():
    pass


def cordex_dataset_id(
    ds,
    sep=".",
    id_attrs=[
        "CORDEX_domain",
        "driving_model_id",
        "institute_id",
        "model_id",
        "experiment_id",
        "frequency",
        #  "driving_model_ensemble_member"
    ],
):
    """Creates a unique string id for e.g. saving files to disk from CORDEX output.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset
    sep : str, optional
        String/Symbol to seperate fields in resulting string, by default "."

    Returns
    -------
    str
        Concatenated

    """
    return _key_from_attrs(ds, id_attrs, sep=sep)
