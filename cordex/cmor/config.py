import datetime as dt

import numpy as np

options = {"table_prefix": "CORDEX-CMIP6", "exit_control": "CMOR_NORMAL"}


# time offsets relative to left labeling for resampling.
# we want the time label to be the center.
loffsets = {
    "3H": dt.timedelta(hours=1, minutes=30),
    "6H": dt.timedelta(hours=3),
    "D": dt.timedelta(hours=12),
}

time_axis_names = {"point": "time1", "mean": "time"}

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


def set_options(**kwargs):
    for k, v in kwargs.items():
        if k in options:
            options[k] = v
        else:
            raise Exception(f"unkown config option: {k}")
