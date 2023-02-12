import tempfile

from . import cf


def get_tempfile():
    """Creates a temporay filename."""
    return tempfile.mkstemp()[1]


def to_center_coordinate(ds):
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    return ds


def pole(ds):
    """Returns rotated pole longitude and latitude"""
    pole_lon = ds[cf.DEFAULT_MAPPING_NCVAR].grid_north_pole_longitude
    pole_lat = ds[cf.DEFAULT_MAPPING_NCVAR].grid_north_pole_latitude
    return pole_lon, pole_lat


def pole_crs(ds):
    """Return a cartropy RotatedPole instance"""
    from cartopy.crs import RotatedPole

    return RotatedPole(*pole(ds))
