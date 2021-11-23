import xarray as xr

from cordex.preprocessing import rename_cordex

from cordex import cordex_domain


def create_test_ds(name, pol_name="rotated_latitude_longitude"):
    domain = cordex_domain(name, mapping_name=pol_name, dummy=True, add_vertices=True)
    return domain


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
    ds = create_wrf_test("EUR-11")
    assert rename_cordex(ds).equals(create_test_ds("EUR-11"))
