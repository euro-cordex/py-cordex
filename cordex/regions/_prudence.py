from warnings import warn

import pandas as pd
import xarray as xr

from ._regions import WGS84
from ._resources import fetch_prudence

message = (
    "Prudence regions are deprecated and will be removed in a future version."
    "Please use regionmask.defined_regions.prudence instead."
)


def _create_polygons(df):
    from shapely.geometry import Polygon

    return [Polygon(_get_vertices(df, area)) for area in df.index]


def _get_vertices(df, area):
    coords = df.loc[area]
    return [
        [coords.west, coords.south],
        [coords.east, coords.south],
        [coords.east, coords.north],
        [coords.west, coords.north],
    ]


def _create_region(df, area):
    import regionmask

    polygon = _get_vertices(df, area)
    return regionmask.Regions([polygon])


def regions(df):
    import regionmask

    regions = [_get_vertices(df, area) for area in df.index]
    return regionmask.Regions(
        regions, names=df.name, abbrevs=df.index, name="prudence regions"
    )


def geodataframe(df, crs=WGS84):
    import geopandas as gpd

    return gpd.GeoDataFrame(
        df.drop(["west", "south", "east", "north"], axis=1),
        index=df.index,
        crs=crs,
        geometry=_create_polygons(df),
    )


def mask_3D(regions, lon, lat, **kwargs):
    regs = regions
    masks = [regs[[i]].mask_3D(lon, lat, **kwargs) for i in regs.numbers]
    return xr.concat(masks, dim="region")


class Prudence:
    """Prudence regions in Europe.

    Prediction of Regional scenarios and Uncertainties for Defining EuropeaN
    Climate change risks and Effects.
    """

    @property
    def df(self):
        warn(
            message,
            DeprecationWarning,
            stacklevel=2,
        )
        return pd.read_csv(fetch_prudence(), na_filter=False, index_col="area")

    @property
    def geodataframe(self):
        warn(
            message,
            DeprecationWarning,
            stacklevel=2,
        )
        return geodataframe(self.df)

    @property
    def regionmask(self):
        warn(
            message,
            DeprecationWarning,
            stacklevel=2,
        )
        return regions(self.df)

    def mask_3D(self, lon, lat, **kwargs):
        warn(
            message,
            DeprecationWarning,
            stacklevel=2,
        )
        return mask_3D(self.regionmask, lon, lat, **kwargs)


prudence = Prudence()
