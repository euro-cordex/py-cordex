

import pandas as pd
import xarray as xr

df = pd.read_csv("https://raw.githubusercontent.com/euro-cordex/tables/master/regions/prudence.csv", na_filter=False, index_col="area")


def _create_polygons(df):
    from shapely.geometry import Polygon
    return [ Polygon(_get_vertices(area)) for area in df.index ]


def _get_vertices(df, area):
    coords = df.loc[area]
    return [[coords.west, coords.south], [coords.east, coords.south], 
            [coords.east, coords.north], [coords.west, coords.north]]


def _create_region(df, area):
    import regionmask
    polygon = _get_vertices(df, area)
    return regionmask.Regions([polygon])


def regions():
    import regionmask
    regions = [ _get_vertices(df, area) for area in df.index ]
    return regionmask.Regions(regions, names=df.name, abbrevs=df.index, name="prudence regions")


def geodataframe(df=df, crs='epsg:4326'):
    import geopadnas as gpd
    return gpd.GeoDataFrame(df.drop(['west', 'south', 'east', 'north'], axis=1), 
                            index=df.index,
                            crs=crs,
                            geometry=_create_polygons(df))


def mask_3D(lon, lat, **kwargs):
    regs = regions()
    masks = [regs[[i]].mask_3D(lon, lat, **kwargs) for i in regs.numbers]
    return xr.concat(masks, dim="region")
