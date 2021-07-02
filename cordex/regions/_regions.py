"""convert to WGS84 latitude-longitude projection"""
WGS84 = "EPSG:4326"


def get_geodataframe(shapefile, to_crs=WGS84, **kwargs):
    import geopandas as gpd

    shp = gpd.read_file(shapefile)
    if to_crs is not None:
        shp = shp.to_crs(to_crs)
    return shp


def get_regionmask(geodataframe, **kwargs):
    import regionmask

    return regionmask.from_geopandas(geodataframe, **kwargs)
