from cartopy import crs as ccrs

crs_map = {
    "rotated_latitude_longitude": {
        "crs": ccrs.RotatedPole,
        "kwargs": {
            "pole_longitude": "grid_north_pole_longitude",
            "pole_latitude": "grid_north_pole_latitude",
            "central_rotated_longitude": 0.0,
        },
    }
}


def _ccrs_from_cf(mapping):
    crs = crs_map.get(mapping.grid_mapping_name)
    kwargs = {kw: (mapping.attrs.get(v) or v) for kw, v in crs["kwargs"].items()}
    return crs["crs"](**kwargs)


def get_ccrs(ds):
    mapping = ds.cf["grid_mapping"]
    return _ccrs_from_cf(mapping)
