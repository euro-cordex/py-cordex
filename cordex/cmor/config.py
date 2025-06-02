import datetime as dt

import numpy as np

options = {
    "table_prefix": "CORDEX-CMIP6",
    "exit_control": "CMOR_NORMAL",
    "resample_kwargs": {"closed": "left"},
    "earth_radius": 6371229.0,  # in meters
}


# time offsets relative to left labeling for resampling.
# we want the time label to be the center.
loffsets = {
    "H": dt.timedelta(minutes=30),
    "3H": dt.timedelta(hours=1, minutes=30),
    "6H": dt.timedelta(hours=3),
    "D": dt.timedelta(hours=12),
}

time_axis_names = {
    "point": "time1",
    "mean": "time",
    "maximum": "time",
    "minimum": "time",
}

units_format = "cf"  # "~P"
# units.define("deg = degree")

# map mip frequencies to pandas frequencies
freq_map = {
    "1hr": "H",
    "1hrPt": "H",
    "3hr": "3H",
    "3hrPt": "3H",
    "6hr": "6H",
    "day": "D",
    "mon": "MS",
}

time_units_default = "days since 1950-01-01T00:00:00"
time_dtype = np.double
time_bounds_name = "time_bounds"
# Y=2000

units_convert_rules = {
    "mm": (lambda x: x * 1.0 / 86400.0, "kg m-2 s-1"),
    "kg/kg": (lambda x: x, "1"),
}

grid_entry_mapping = {
    "rotated_latitude_longitude": {
        "X": "grid_longitude",
        "Y": "grid_latitude",
        "default_attrs": {
            "grid_north_pole_latitude": [None, ""],
            "grid_north_pole_longitude": [None, ""],
            "north_pole_grid_longitude": [0.0, ""],
            "earth_radius": [options["earth_radius"], ""],
        },
    },
    "lambert_conformal_conic": {
        "X": "x",
        "Y": "y",
        "default_attrs": {
            "standard_parallel": [None, ""],
            "longitude_of_central_meridian": [None, ""],
            "latitude_of_projection_origin": [None, ""],
            "false_easting": [0.0, ""],
            "false_northing": [0.0, ""],
            "earth_radius": [options["earth_radius"], ""],
        },
    },
}


def set_options(**kwargs):
    for k, v in kwargs.items():
        if k in options:
            options[k] = v
        else:
            raise Exception(f"unkown config option: {k}")
