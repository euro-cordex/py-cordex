from cartopy import crs as ccrs
from .utils import get_grid_mapping_attrs
from pyproj import CRS
from xarray import Dataset

crs_map = {
    "rotated_latitude_longitude": {
        "crs": ccrs.RotatedPole,
        "kwargs": {
            "pole_longitude": "grid_north_pole_longitude",
            "pole_latitude": "grid_north_pole_latitude",
            "central_rotated_longitude": 0.0,
        },
    }
}


def _ccrs_from_cf(attrs):
    crs = crs_map.get(attrs["grid_mapping_name"])
    kwargs = {kw: (attrs.get(v) or v) for kw, v in crs["kwargs"].items()}
    return crs["crs"](**kwargs)


def get_ccrs(obj):
    if isinstance(obj, Dataset):
        attrs = obj.cf["grid_mapping"].attrs
    elif isinstance(obj, str):
        attrs = get_grid_mapping_attrs(obj)
    return _ccrs_from_cf(attrs)


def get_crs(domain_id):
    """creates a pyproj CRS instance from a CORDEX domain_id"""
    # i go via from_cf to get the correct crs
    attrs = get_grid_mapping_attrs(domain_id)
    return CRS.from_cf(attrs)
    # proj4 = f"+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p={dm.pollat} +lon_0={180+dm.pollon} +datum=WGS84 +no_defs +type=crs"
