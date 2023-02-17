import numpy as np
import pytest
import xarray as xr

import cordex as cx
from cordex import cordex_domain
from cordex.preprocessing.preprocessing import (
    attr_to_coord,
    check_domain,
    cordex_dataset_id,
    get_grid_mapping,
    get_grid_mapping_name,
    member_id_to_dset_id,
    promote_empty_dims,
    remap_lambert_conformal,
    rename_cordex,
    replace_coords,
    sort_ds_dict_by_attr,
)

from . import requires_xesmf


def create_test_ds(
    name, pol_name="rotated_latitude_longitude", dummy=True, bounds=True, **kwargs
):
    domain = cordex_domain(
        name, mapping_name=pol_name, dummy=dummy, bounds=bounds, **kwargs
    )
    return domain


@pytest.fixture
def test_ensemble():
    return cx.tutorial.ensemble()


def create_wrf_test(name):
    ds = create_test_ds("EUR-11", "rotated_pole")
    # reindex dummy with lon lat like in wrf
    ds["dummy"] = xr.DataArray(
        ds.dummy.values, dims=("lat", "lon"), attrs=ds.dummy.attrs
    )
    return ds


# def create_test_ds(xname, yname, xlen, ylen, name):
#     x = np.linspace(0, 359, xlen)
#     y = np.linspace(-90, 89, ylen)

#     data = np.random.rand(len(x), len(y))
#     ds = xr.DataArray(data, coords=[(xname, x), (yname, y)]).to_dataset(
#         name="test"
#     )
#     ds.attrs["source_id"] = "test_id"
#     # if x and y are not lon and lat, add lon and lat to make sure there are no conflicts
#     lon = ds[xname] * xr.ones_like(ds[yname])
#     lat = xr.ones_like(ds[xname]) * ds[yname]
#     if xname != "lon" and yname != "lat":
#         ds = ds.assign_coords(lon=lon, lat=lat)
#     else:
#         ds = ds.assign_coords(longitude=lon, latitude=lat)
#     return ds


def test_wrf_case():
    """Test the wrf exception"""
    ds = xr.decode_cf(rename_cordex(create_wrf_test("EUR-11")), decode_coords="all")
    xr.testing.assert_equal(
        ds, xr.decode_cf(create_test_ds("EUR-11"), decode_coords="all")
    )


@pytest.mark.parametrize("lon_name", ["longitude"])
@pytest.mark.parametrize("lat_name", ["latitude"])
@pytest.mark.parametrize("pol_name", ["rotated_latitude_longitude", "rotated_pole"])
@pytest.mark.parametrize("lon_vertices", ["longitude_vertices"])
@pytest.mark.parametrize("lat_vertices", ["latitude_vertices"])
def test_rename_cordex(lon_name, lat_name, pol_name, lon_vertices, lat_vertices):
    dm = create_test_ds("EUR-11", pol_name)
    dm = dm.rename(
        {
            "lon": lon_name,
            "lat": lat_name,
            "lon_vertices": lon_vertices,
            "lat_vertices": lat_vertices,
        }
    )
    # assert rename_cordex(dm).equals(create_test_ds("EUR-11"))
    xr.testing.assert_equal(
        xr.decode_cf(rename_cordex(dm), decode_coords="all"),
        xr.decode_cf(create_test_ds("EUR-11"), decode_coords="all"),
    )


def test_grid_mapping():
    ds = create_test_ds("EUR-11")
    assert get_grid_mapping(ds).equals(ds.rotated_latitude_longitude)


def test_replace_coords():
    ds = create_test_ds("EUR-11")
    ds["rlon"] = np.arange(ds.rlon.size)
    ds["rlat"] = np.arange(ds.rlat.size)
    ds["lon"] = np.arange(ds.lon.size)
    ds["lat"] = np.arange(ds.lon.size)
    assert replace_coords(ds).equals(create_test_ds("EUR-11"))


def test_cordex_dataset_id():
    ds = create_test_ds("EUR-11", attrs="CORDEX")
    ds.attrs["driving_model_id"] = "MY-DRIVE-MODEL"
    ds.attrs["institute_id"] = "INSTITUTE"
    ds.attrs["model_id"] = "RCM"
    ds.attrs["experiment_id"] = "historical"
    ds.attrs["frequency"] = "mon"
    assert (
        cordex_dataset_id(ds, sep=".")
        == "EUR-11.MY-DRIVE-MODEL.INSTITUTE.RCM.historical.mon"
    )


def test_member_id_to_dset_id():
    ds = create_test_ds("EUR-11", attrs="CORDEX")
    ds.attrs["driving_model_id"] = "MY-DRIVE-MODEL"
    ds.attrs["institute_id"] = "INSTITUTE"
    ds.attrs["model_id"] = "RCM"
    ds.attrs["experiment_id"] = "historical"
    ds.attrs["frequency"] = "mon"
    ds.attrs["member"] = "r1i1p1"

    ds = attr_to_coord(ds, "member")
    ds_id_old = cordex_dataset_id(
        ds, sep="."
    )  # == 'EUR-11.MY-DRIVE-MODEL.INSTITUTE.RCM.historical.mon'
    ds_dict = {ds_id_old: ds}
    ds_dict_new = member_id_to_dset_id(ds_dict)
    for ds_id, ds in ds_dict_new.items():
        assert ds_id == ds_id_old + ".{}".format(ds.attrs["member"])


def test_promote_empty_dims():
    ds = create_test_ds("EUR-11")
    ds = ds.drop_vars(["rlon", "rlat"])
    ds_promoted = promote_empty_dims(ds)
    assert set(["rlon", "rlat"]).issubset(set(ds_promoted.coords))


def test_sort_ds_dict_by_attr(test_ensemble):
    ds_dict_sorted = sort_ds_dict_by_attr(test_ensemble, "experiment_id")
    assert "evaluation" in ds_dict_sorted


def test_check_domain():
    assert check_domain(cordex_domain("EUR-11"))


@requires_xesmf
def test_remap_lambert_conformal(test_ensemble):
    remap = {key: remap_lambert_conformal(ds) for key, ds in test_ensemble.items()}
    for ds in remap.values():
        assert get_grid_mapping_name(ds) in [
            "rotated_latitude_longitude",
            "rotated_pole",
        ]
