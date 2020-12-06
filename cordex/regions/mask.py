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




class GERMANY():

    @staticmethod
    def VG2500(domain='lan'):
        url = addr.VG2500(domain)
        geodata = get_geodata(url)
        geodata['mask_name'] = geodata["ARS"] + '_' + geodata['GEN']
        return get_regionmask(geodata, names="mask_name", abbrevs='_from_name')

