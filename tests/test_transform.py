import numpy as np
import xarray as xr
from cf_xarray.datasets import rotds
from pyproj import CRS

from cordex import (
    derotate_vector,
    domain,
    rotated_coord_transform,
    map_crs,
    transform,
    vertices,
    transform_bounds,
    transform_coords,
)

from . import requires_cartopy


@requires_cartopy
def test_map_crs():
    import cartopy.crs as ccrs

    eur11 = domain("EUR-11")
    pole = (
        eur11.rotated_latitude_longitude.grid_north_pole_longitude,
        eur11.rotated_latitude_longitude.grid_north_pole_latitude,
    )
    lon1, lat1 = rotated_coord_transform(
        eur11.rlon, eur11.rlat, *pole, direction="rot2geo"
    )
    transform = ccrs.RotatedPole(*pole)
    lon2, lat2 = map_crs(eur11.rlon, eur11.rlat, transform)

    assert np.allclose(lon1.T, lon2)
    assert np.allclose(lat1.T, lat2)

    # test if retransforming of lon lat to rlon rlat gives correct results
    rlon2, rlat2 = map_crs(
        eur11.lon, eur11.lat, src_crs=ccrs.PlateCarree(), trg_crs=transform
    )
    rlat1, rlon1 = xr.broadcast(eur11.rlat, eur11.rlon)

    assert np.allclose(rlon1, rlon2)
    assert np.allclose(rlat1, rlat2)


def test_transform():
    ds = rotds.copy(deep=False)
    rotated = CRS.from_cf(ds.cf["grid_mapping"].attrs)

    xt, yt = transform(ds.rlon, ds.rlat, rotated)
    assert np.allclose(ds.lon, xt)
    assert np.allclose(ds.lat, yt)

    # test if retransforming of lon lat to rlon rlat gives correct results
    rlon_t, rlat_t = transform(xt, yt, CRS("EPSG:4326"), rotated)
    rlat_b, rlon_b = xr.broadcast(ds.rlat, ds.rlon)

    assert np.allclose(rlon_b, rlon_t)
    assert np.allclose(rlat_b, rlat_t)

    ds1 = transform_coords(ds)

    assert np.allclose(ds1.lon, ds1.xt)
    assert np.allclose(ds1.lat, ds1.yt)


def test_derotate_vector():
    u = np.array([1.0, 0.0, -1.0])
    v = np.array([-1.0, 1.0, 0])
    lon = np.array([0.0, 0.0, 0.0])
    lat = np.array([0.0, 0.0, 0.0])

    pollon = 90.0
    pollat = 45.0
    u1_expect = np.array([0.0, 0.70710678, -0.70710678])
    v1_expect = np.array([-1.41421356, 0.70710678, 0.70710678])
    u1, v1 = derotate_vector(u, v, lon, lat, pollon, pollat)
    assert np.allclose(u1, u1_expect)
    assert np.allclose(v1, v1_expect)

    pollon = 180.0
    pollat = 90.0
    u1_expect = u
    v1_expect = v
    u1, v1 = derotate_vector(u, v, lon, lat, pollon, pollat)
    assert np.allclose(u1, u1_expect)
    assert np.allclose(v1, v1_expect)


def test_bounds():
    # assert that we get the same bounds as before
    ds = domain("EUR-11", bounds=True)
    _v = vertices(ds.rlon, ds.rlat, CRS.from_cf(ds.cf["grid_mapping"].attrs))
    v = transform_bounds(ds)
    np.array_equal(v.lon_vertices, _v.lon_vertices)
    np.array_equal(v.lat_vertices, _v.lat_vertices)

    assert "longitude" in v.cf.bounds
    assert "latitude" in v.cf.bounds
