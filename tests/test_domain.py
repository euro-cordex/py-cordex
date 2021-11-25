# -*- coding: utf-8 -*-
# flake8: noqa

import numpy as np
import pytest

import cordex as cx

from . import has_cartopy, requires_cartopy


def test_domain_basic():
    eur11 = cx.cordex_domain("EUR-11")
    assert "lon_vertices" in cx.cordex_domain("EUR-11", add_vertices=True)
    assert "lat_vertices" in cx.cordex_domain("EUR-11", add_vertices=True)
    assert "rotated_pole" in cx.cordex_domain("EUR-11", mapping_name="rotated_pole")
    assert "rotated_pole" in cx.cordex_domain("EUR-11", mapping_name="rotated_pole")
    assert "dummy" in cx.cordex_domain("EUR-11", dummy=True)
    assert "topo" in cx.cordex_domain("EUR-11", dummy="topo")
    assert cx.cordex_domain("EUR-11").attrs["CORDEX_domain"] == "EUR-11"
    assert "institution" in cx.cordex_domain("EUR-11", attrs="CORDEX").attrs


def test_constructor():
    eur11 = cx.cordex_domain("EUR-11")
    eur11_user = cx.create_dataset(
        nlon=424,
        nlat=412,
        dlon=0.11,
        dlat=0.11,
        ll_lon=-28.375,
        ll_lat=-23.375,
        pollon=-162.00,
        pollat=39.25,
    )
    assert eur11_user.equals(eur11)


def test_domain_info():
    import pandas as pd

    info = {
        "short_name": "EUR-11",
        "region": 4,
        "long_name": "Europe",
        "nlon": 424,
        "nlat": 412,
        "ll_lon": -28.375,
        "ur_lon": 18.155,
        "ll_lat": -23.375,
        "ur_lat": 21.835,
        "dlon": 0.11,
        "dlat": 0.11,
        "pollon": -162.0,
        "pollat": 39.25,
    }
    table = pd.DataFrame.from_dict(
        {key: [item] for key, item in info.items()}
    ).set_index("short_name")
    assert cx.domain_info("EUR-11", table) == info


# @requires_cartopy
def test_mapping():
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

    assert np.allclose(lon1, lon2)
    assert np.allclose(lat1, lat2)


def test_vertices():
    eur11 = cx.cordex_domain("EUR-11")
    import cartopy.crs as ccrs

    pole = (
        eur11.rotated_latitude_longitude.grid_north_pole_longitude,
        eur11.rotated_latitude_longitude.grid_north_pole_latitude,
    )
    vertices = cx.vertices(eur11.rlon, eur11.rlat, src_crs=ccrs.RotatedPole(*pole))
