import numpy as np
import xarray as xr
from cf_xarray.datasets import rotds
from pyproj import CRS

import cordex as cx

from . import requires_cartopy


@requires_cartopy
def test_map_crs():
    import cartopy.crs as ccrs

    eur11 = cx.cordex_domain("EUR-11")
    pole = (
        eur11.rotated_latitude_longitude.grid_north_pole_longitude,
        eur11.rotated_latitude_longitude.grid_north_pole_latitude,
    )
    lon1, lat1 = cx.rotated_coord_transform(
        eur11.rlon, eur11.rlat, *pole, direction="rot2geo"
    )
    transform = ccrs.RotatedPole(*pole)
    lon2, lat2 = cx.map_crs(eur11.rlon, eur11.rlat, transform)

    assert np.allclose(lon1.T, lon2)
    assert np.allclose(lat1.T, lat2)

    # test if retransforming of lon lat to rlon rlat gives correct results
    rlon2, rlat2 = cx.map_crs(
        eur11.lon, eur11.lat, src_crs=ccrs.PlateCarree(), trg_crs=transform
    )
    rlat1, rlon1 = xr.broadcast(eur11.rlat, eur11.rlon)

    assert np.allclose(rlon1, rlon2)
    assert np.allclose(rlat1, rlat2)


def test_transform():
    ds = rotds.copy(deep=False)
    rotated = CRS.from_cf(ds.cf["grid_mapping"].attrs)

    xt, yt = cx.transform(ds.rlon, ds.rlat, rotated)
    assert np.allclose(ds.lon, xt)
    assert np.allclose(ds.lat, yt)

    # test if retransforming of lon lat to rlon rlat gives correct results
    rlon_t, rlat_t = cx.transform(xt, yt, CRS("EPSG:4326"), rotated)
    rlat_b, rlon_b = xr.broadcast(ds.rlat, ds.rlon)

    assert np.allclose(rlon_b, rlon_t)
    assert np.allclose(rlat_b, rlat_t)
