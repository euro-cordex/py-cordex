"""This module should provide some help to mask data.
"""

import geopandas as gpd
import regionmask

from . import _address as addr

"""convert to WGS84 latitude-longitude projection"""
WGS84 = "EPSG:4326"


def get_geodata(shapefile, to_crs=WGS84, **kwargs):
    shp = gpd.read_file(shapefile)
    if to_crs is not None:
        shp = shp.to_crs(to_crs)
    return shp


def get_regionmask(geodata, **kwargs):
    return regionmask.from_geopandas(geodata, **kwargs)


# tas_mask = laender.mask(tas, lon_name='lon', lat_name='lat')
# mask_2D = regionmask.mask_geopandas(shp, tas.lon, tas.lat)


def gridded_mask(rmask, da=None, lon=None, lat=None):
    if da is not None:
        if lon is None:
            lon = "lon"
        if lat is None:
            lat = "lat"
        return rmask.mask(da, lon_name=lon, lat_name=lat)
    else:
        if lon is None or lat is None:
            raise Exception(
                "Missing lon or lat attribute for regionmask.mask_geopandas"
            )
        return regionmask.mask_geopandas(rmask, lon, lat)
