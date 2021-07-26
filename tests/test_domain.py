# -*- coding: utf-8 -*-
# flake8: noqa

import numpy as np
import pytest

import cordex as cx

from . import has_cartopy, requires_cartopy


def test_constructor():
    eur11 = cx.cordex_domain('EUR-11')
    eur11_user = cx.create_dataset(nlon=424, nlat=412, dlon=0.11, dlat=0.11, ll_lon=-28.375, ll_lat=-23.375, pollon=-162.00, pollat=39.25)
    assert(eur11_user.equals(eur11))


@requires_cartopy
def test_mapping():
    import cartopy.crs as ccrs
    eur11 = cx.cordex_domain('EUR-11')
    pole = ( eur11.rotated_latitude_longitude.grid_north_pole_longitude, eur11.rotated_latitude_longitude.grid_north_pole_latitude )
    lon1, lat1 = cx.rotated_coord_transform(eur11.rlon, eur11.rlat, *pole, direction="rot2geo")
    transform = ccrs.RotatedPole(-162.0, 39.25)
    lon2, lat2 = cx.map_crs(eur11.rlon, eur11.rlat, transform)

    assert(np.allclose(lon1, lon2))
    assert(np.allclose(lat1, lat2))

