import os
from os import path as op
from warnings import warn
from collections.abc import Iterable
from collections import OrderedDict

import cf_xarray as cfxr
import pandas as pd
import xarray as xr

# from cf_xarray.units import units

try:
    import cmor
except Exception:
    warn("no python cmor package available, consider installing it")

# import cordex as cx
from .. import create_dataset
from ..domain import domain
from .config import (
    freq_map,
    loffsets,
    options,
    time_axis_names,
    time_dtype,
    time_units_default,
    units_format,
    grid_entry_mapping,
)

# from .derived import derivator
from .utils import (
    _encode_time,
    _get_cfvarinfo,
    _get_cordex_pole,
    _read_table,
    _strip_time_cell_method,
    _tmp_table,
    mid_of_month,
    month_bounds,
    time_bounds_name,
)

xr.set_options(keep_attrs=True)

try:
    import flox  # noqa

    has_flox = True
    default_flox_method = "blockwise"
except ImportError:
    has_flox = False
    default_flox_method = None


def resample_both_closed(ds, hfreq, op, **kwargs):
    rolling = getattr(ds.rolling(time=hfreq + 1, center=True), op)()
    freq = f"{hfreq}H"
    return rolling.resample(time=freq, loffset=0.5 * pd.Timedelta(hfreq, "H")).nearest()


def _resample_op(ds, hfreq, op, **kwargs):
    rolling = getattr(ds.rolling(time=hfreq + 1, center=True), op)()
    freq = f"{hfreq}H"
    return rolling.resample(time=freq, loffset=0.5 * pd.Timedelta(hfreq, "H")).nearest()


def _get_loffset(time):
    return loffsets.get(time, None)


def _clear_time_axis(ds):
    """Delete timesteps with NaN arrays"""
    for data_var in ds.data_vars:
        ds = ds.dropna(dim="time", how="all")
    return ds


def _resample(
    ds, time, time_cell_method="point", label="left", time_offset=True, **kwargs
):
    """Resample a variable."""
    if time_cell_method == "point":
        return ds.resample(
            time=time, label="right", closed="right"
        ).nearest()  # .interpolate("nearest") # use as_freq?
    elif time_cell_method == "mean":
        if time_offset is True:
            loffset = _get_loffset(time)
        else:
            loffset = None
        mean_kwargs = {}
        if has_flox:
            mean_kwargs["engine"] = "flox"
            mean_kwargs["method"] = default_flox_method
        return ds.resample(time=time, label=label, offset=loffset, **kwargs).mean(
            **mean_kwargs
        )
    else:
        raise Exception(f"unknown time_cell_method: {time_cell_method}")


def _get_bnds(values):
    bnds = [None] * (len(values) + 1)
    bnds[0] = values[0] - (values[1] - values[0]) / 2
    bnds[len(values)] = values[-1] + (values[-1] - values[-2]) / 2
    i = 1
    while i < len(values):
        bnds[i] = values[i] - (values[i] - values[i - 1]) / 2
        i += 1
    return bnds


def _crop_to_cordex_domain(ds, domain_id, tolerance=1.0e-2):
    """Crop data to official CORDEX domain"""
    try:
        grid = domain(domain_id)
    except KeyError:
        raise KeyError(
            f'unknown domain_id: "{domain_id}", if you are using an unofficial domain, set crop=False'
        )

    info = grid.cx.info()
    rlon_tol = tolerance * info["dlon"]
    rlat_tol = tolerance * info["dlat"]
    # the method=='nearest' approach does not work well with dask
    return ds.sel(
        rlon=slice(grid.rlon.min() - rlon_tol, grid.rlon.max() + rlon_tol),
        rlat=slice(grid.rlat.min() - rlat_tol, grid.rlat.max() + rlat_tol),
    )


def _load_table(table):
    cmor.load_table(table)


def _setup(dataset_table, mip_table, grids_table=None, inpath="."):
    if grids_table is None:
        grids_table = f'{options["table_prefix"]}_grids.json'
    cmor.setup(
        inpath,
        set_verbosity=cmor.CMOR_NORMAL,
        netcdf_file_action=cmor.CMOR_REPLACE,
        exit_control=getattr(cmor, options["exit_control"]),
        logfile=None,
    )
    cmor.dataset_json(dataset_table)
    grid_id = cmor.load_table(grids_table)
    table_id = cmor.load_table(mip_table)
    cmor.set_table(table_id)
    return {"grid": grid_id, "mip": table_id}


def _get_time_axis_name(time_cell_method):
    """Get the name of the CMOR time coordinate"""
    return time_axis_names.get(time_cell_method, "time")


def _define_grid(ds, table_id, rewrite_grid="auto"):
    if rewrite_grid == "auto":
        try:
            grid_attrs = ds.cx.info()
            grid = create_dataset(**grid_attrs, bounds=True)
        except (KeyError, ValueError):
            warn("can not rewrite grid")
            grid = ds
    else:
        grid = ds
    grid_mapping = grid.cf["grid_mapping"]
    grid_mapping_name = grid_mapping.grid_mapping_name
    entry_mapping = grid_entry_mapping.get(grid_mapping_name)
    default_attrs = entry_mapping["default_attrs"]
    # if "domain_id" in ds.attrs:
    #     try:
    #         grid = domain(ds.attrs["domain_id"], bounds=True)
    #         lon_vertices = grid.lon_vertices.to_numpy()
    #         lat_vertices = grid.lat_vertices.to_numpy()
    #     except KeyError:
    #         warn(f"Unknown domain: {ds.attrs['domain_id']}")
    # else:
    #     lon_vertices = None
    #     lat_vertices = None

    # if "longitude" not in ds.cf.coords or "latitude" not in ds.cf.coords:
    #     ds = cx.transform_coords(ds, trg_dims=("lon", "lat"))

    cmor.set_table(table_id)

    cmorY = cmor.axis(
        table_entry=entry_mapping["Y"],
        coord_vals=grid.cf["Y"].to_numpy(),
        units=grid.cf["Y"].units,
    )
    cmorX = cmor.axis(
        table_entry=entry_mapping["X"],
        coord_vals=grid.cf["X"].to_numpy(),
        units=grid.cf["X"].units,
    )

    latitude = grid.cf["latitude"].to_numpy()
    longitude = grid.cf["longitude"].to_numpy()
    latitude_vertices = (
        grid.cf.get_bounds("latitude").to_numpy()
        if "latitude" in grid.cf.bounds
        else None
    )
    longitude_vertices = (
        grid.cf.get_bounds("longitude").to_numpy()
        if "longitude" in grid.cf.bounds
        else None
    )

    cmorGrid = cmor.grid(
        [cmorY, cmorX],
        latitude=latitude,
        longitude=longitude,
        latitude_vertices=latitude_vertices,
        longitude_vertices=longitude_vertices,
    )

    # attrs_dict = {k: grid_mapping.attrs.get(v) or v for k, v in default_attrs.items()}
    attrs_dict = OrderedDict()
    for k, v in default_attrs.items():
        if k in grid_mapping.attrs:
            value = grid_mapping.attrs[k]
            if isinstance(value, Iterable):
                for i, val in enumerate(value):
                    # see
                    attrs_dict[k + str(i + 1)] = [val, v[1]]
                    # attrs_dict[k] = [val, v[1]]
            else:
                attrs_dict[k] = [value, v[1]]
        elif v[0] is not None:
            attrs_dict[k] = v
        else:
            raise KeyError(
                f"Missing attribute {k} in grid mapping variable: {grid_mapping_name}"
            )

    cmor.set_grid_mapping(
        cmorGrid,
        grid_mapping_name,
        # parameter_names=attrs_dict,
        parameter_names=list(attrs_dict.keys()),
        parameter_values=[v[0] for v in attrs_dict.values()],
        parameter_units=[v[1] for v in attrs_dict.values()],
        # len(attrs_dict) * [""],
    )

    return cmorGrid


def _time_bounds(ds, freq=None):
    return cfxr.bounds_to_vertices(ds.cf.get_bounds("time"), "bounds")


def _define_time(ds, table_id, time_cell_method=None):
    cmor.set_table(table_id)

    if time_cell_method is None:
        warn("no time_cell_method given, assuming: mean")
        time_cell_method = "mean"

    # encode time and time bounds

    time_axis_encode = _encode_time(ds.time).to_numpy()
    time_axis_name = _get_time_axis_name(time_cell_method)

    if time_cell_method != "point":
        time_bounds = _time_bounds(ds)
        time_bounds.encoding = ds.time.encoding
        time_bounds_encode = _encode_time(time_bounds).to_numpy()
    else:
        time_bounds_encode = None

    return cmor.axis(
        time_axis_name,
        coord_vals=time_axis_encode,
        # cell_bounds=_get_bnds(time_encode.to_numpy()),
        cell_bounds=time_bounds_encode,
        units=ds.time.encoding["units"],
    )


def _define_axis(ds, table_ids, time_cell_method="point"):
    cmorGrid = _define_grid(ds, table_ids["grid"])

    if "time" in ds:
        cmorTime = _define_time(ds, table_ids["mip"], time_cell_method)
    else:
        cmorTime = None

    # add z axis if required
    if "Z" in ds.cf.dims:
        z = ds.cf["Z"]
        print(ds.cf.bounds)
        if "Z" not in ds.cf.bounds:
            # try to add bounds via cf_xarray
            warn("found not bounds for z-axis, trying to create bounds...")
            ds = ds.cf.add_bounds("Z", output_dim="bnds")
        bounds = ds.cf.get_bounds("Z")
        if bounds.ndim != 1:
            # transpose to bounds
            bounds = cfxr.bounds_to_vertices(bounds, "bnds")
        cmorZ = cmor.axis(
            table_entry=z.name,
            coord_vals=z.to_numpy(),
            units=z.attrs.get("units"),
            cell_bounds=bounds.to_numpy(),
        )
    else:
        cmorZ = None

    return cmorTime, cmorZ, cmorGrid


def _cmor_write(da, table_id, cmorTime, cmorZ, cmorGrid, file_name=True):
    """write to netcdf via cmor python API"""

    cmor.set_table(table_id)

    # create coordinate ids
    coords = []
    if cmorTime:
        coords.append(cmorTime)
    if cmorZ:
        coords.append(cmorZ)
    coords.append(cmorGrid)

    cmor_var_kwargs = {}
    for kwarg in ["positive", "missing_value", "original_name", "history", "comment"]:
        if kwarg in da.attrs:
            cmor_var_kwargs[kwarg] = da.attrs[kwarg]

    cmor_var = cmor.variable(
        table_entry=da.name, units=da.units, axis_ids=coords, **cmor_var_kwargs
    )

    if "time" in da.coords:
        ntimes_passed = da.time.size
    else:
        ntimes_passed = None
    cmor.write(cmor_var, da.to_numpy(), ntimes_passed=ntimes_passed)

    return cmor.close(cmor_var, file_name=file_name)


def _units_convert(da, cf_units, format=None):
    import pint_xarray  # noqa
    from cf_xarray.units import units  # noqa

    if format is None:
        format = units_format
    da_quant = da.pint.quantify()
    da = da_quant.pint.to(cf_units).pint.dequantify(format=units_format)
    return da


def _cf_units_convert(da, table, mapping_table={}):
    """Convert units.

    Convert units according to the rules in units_convert_rules dict.
    Maybe metpy can do this also: https://unidata.github.io/MetPy/latest/tutorials/unit_tutorial.html

    """
    from cf_xarray.units import units as cfxr_units  # noqa

    if da.name in mapping_table:
        map_units = mapping_table[da.name].get("units")
        atr_units = da.attrs.get("units")
        if map_units is not None and atr_units is not None and atr_units != map_units:
            warn(
                f"unit [{map_units}] from mapping table differs from units attribute [{da}], assuming [{map_units}] is correct"
            )
            units = map_units
        elif map_units is not None:
            units = map_units
        else:
            units = atr_units
    else:
        try:
            units = da.units
        except AttributeError:
            warn(
                f"could not find units attribute on {da.name}, if you are sure about the units, set allow_units_convert=False."
            )
    da.attrs["units"] = units
    cf_units = table["variable_entry"][da.name]["units"]

    if cfxr_units.Unit(units) != cfxr_units.Unit(cf_units):
        warn(f"converting units {units} from input data to CF units {cf_units}")
        da = _units_convert(da, cf_units)
        # da = da.pint.to(cf_units)
    return da


# def _convert_cmor_to_resample_frequency(cmor_table):
#    """Convert CMOR table name into resample frequency"""
#    return resample_frequency[cmor_table]


def _get_time_units(ds):
    """Determine time units of dataset"""
    try:
        return ds.time.encoding["units"]
    except Exception:
        return ds.time.units
    return None


def _set_time_encoding(ds, units, orig):
    u = None
    if units is not None:
        if units == "input":
            if "units" in orig.time.encoding:
                u = orig.time.encoding["units"]
            elif "units" in orig.time.attrs:
                u = orig.time.attrs["units"]
        else:
            u = units
    if u is None:
        u = time_units_default
        warn(f"time units are set to default: {u}")
    if not isinstance(u, str):
        raise TypeError(f"time units invalid: {u}")
    ds.time.encoding["units"] = u
    ds.time.encoding["dtype"] = time_dtype
    return ds


def _update_time_axis(ds, freq=None, time_cell_method=None):
    if time_cell_method is None:
        warn("no time_cell_method given, assuming: point")
        time_cell_method = "point"

    if time_cell_method == "mean":
        time_bounds = _time_bounds(ds)
        time_bounds.encoding = ds.time.encoding
        time_bounds_encode = _encode_time(time_bounds).to_numpy()
    else:
        time_bounds_encode = None
    return time_bounds_encode


def _rewrite_time_axis(ds, freq=None, calendar=None):
    """rewrite time axis to ensure correct timestamps

    For monthly frequencies, the timestamp will be mid of month.

    """
    ds = ds.copy(deep=False)
    pd_freq = freq_map.get(freq, None)
    if freq is None:
        pd_freq = xr.infer_freq(ds.time).upper()
    if calendar is None:
        calendar = ds.time.dt.calendar
    start = ds.time.data[0]
    if freq == "mon":
        ds["time"] = mid_of_month(ds)
        return ds

    date_range = xr.cftime_range(
        start,
        periods=ds.time.size,
        freq=pd_freq,
        calendar=calendar,  # inclusive="left"
    )

    return ds.assign_coords(time=date_range)


def _add_month_bounds(ds):
    ds[time_bounds_name] = month_bounds(ds)
    return ds


def _add_time_bounds(ds, cf_freq):
    """add time bounds

    Take special care of monthly frequencies.

    """
    # monthly time bounds are funny in ESGF, it seems that
    # they should always be the first of each month and first
    # of second month. This is not really the bounds you would get
    # from arithmetics but seems fine. We have to take special care
    # of them here.
    if cf_freq == "mon":
        ds = _add_month_bounds(ds)
    else:
        try:
            ds = ds.convert_calendar(
                ds.time.dt.calendar, use_cftime=False
            ).cf.add_bounds("time")
        except Exception:
            # wait for cftime arithemtics in xarry here:
            warn("could not add time bounds.")

    ds[time_bounds_name].encoding = ds.time.encoding
    ds.time.attrs.update({"bounds": time_bounds_name})
    return ds


def _adjust_frequency(ds, cf_freq, input_freq=None, time_cell_method=None):
    if input_freq is None and "time" in ds.coords:
        input_freq = xr.infer_freq(ds.time)
    if input_freq is None:
        warn("Could not determine frequency of input data, will assume it is correct.")
        return ds
    else:
        input_freq = input_freq.upper()
    pd_freq = freq_map[cf_freq]
    if pd_freq != input_freq:
        warn(f"resampling input data from {input_freq} to {pd_freq}")
        resample = _resample(
            ds, pd_freq, time_cell_method=time_cell_method, **options["resample_kwargs"]
        )
        return resample
    return ds


def cmorize_cmor(
    ds, out_name, cmor_table, dataset_table, grids_table=None, inpath=None
):
    # get meta info from cmor table
    if inpath is None:
        inpath = os.path.dirname(cmor_table)

    dataset_table_json = dataset_table
    cmor_table_json = cmor_table

    if isinstance(dataset_table, dict):
        dataset_table_json = _tmp_table(dataset_table)
    if isinstance(cmor_table, dict):
        cmor_table_json = _tmp_table(cmor_table)

    cfvarinfo = _get_cfvarinfo(out_name, cmor_table)

    if cfvarinfo is None:
        raise Exception(f"{out_name} not found in {cmor_table}")

    time_cell_method = _strip_time_cell_method(cfvarinfo)

    table_ids = _setup(
        dataset_table_json, cmor_table_json, grids_table=grids_table, inpath=inpath
    )

    cmorTime, cmorZ, cmorGrid = _define_axis(
        ds, table_ids, time_cell_method=time_cell_method
    )

    return _cmor_write(ds[out_name], table_ids["mip"], cmorTime, cmorZ, cmorGrid)


def prepare_variable(
    ds,
    out_name,
    cmor_table,
    mapping_table=None,
    replace_coords=False,
    allow_units_convert=False,
    allow_resample=False,
    input_freq=None,
    domain_id=None,
    time_units=None,
    rewrite_time_axis=False,
    use_cftime=False,
    squeeze=True,
    crop=None,
    guess_coord_axis=None,
):
    """prepares a variable for cmorization."""

    if mapping_table is None:
        mapping_table = {}

    ds = ds.copy(deep=False)
    # use cf_xarray to guess coordinate meta data
    if guess_coord_axis is None:
        guess_coord_axis = "X" not in ds.cf.dims or "Y" not in ds.cf.dims
    if guess_coord_axis is True:
        ds = ds.cf.guess_coord_axis(verbose=True)

    if isinstance(cmor_table, str):
        cmor_table = _read_table(cmor_table)
    cfvarinfo = _get_cfvarinfo(out_name, cmor_table)

    cf_freq = cfvarinfo["frequency"]
    time_cell_method = _strip_time_cell_method(cfvarinfo)

    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()

    # ensure that we propagate everything
    # ds = xr.decode_cf(ds, decode_coords="all")

    # no mapping table provided, we assume datasets has already correct out_names and units.
    if out_name in ds.data_vars:
        var_ds = ds  # [[out_name]]
    elif mapping_table and out_name not in mapping_table:
        raise Exception(
            f"Could not find {out_name} in dataset. Please make sure, variable names and units have CF standard or pass a mapping table."
        )
    else:
        varname = mapping_table[out_name]["varname"]
        # cf_name = varinfo["cf_name"]
        var_ds = ds  # [[varname]]  # .to_dataset()
        var_ds = var_ds.rename({varname: out_name})
    # remove point coordinates, e.g, height2m
    if squeeze is True:
        var_ds = var_ds.squeeze(drop=True)
    if crop is True:
        # var_ds.attrs["domain_id"] = domain_id
        var_ds = _crop_to_cordex_domain(var_ds, domain_id)
    if replace_coords is True:
        # domain_id = domain_id or var_ds.cx.domain_id
        grid = domain(domain_id, bounds=True)
        var_ds = var_ds.assign_coords(rlon=grid.rlon, rlat=grid.rlat)
        var_ds = var_ds.assign_coords(lon=grid.lon, lat=grid.lat)
        var_ds = var_ds.assign_coords(
            lon_vertices=grid.lon_vertices, lat_vertices=grid.lat_vertices
        )

    if "time" in var_ds:
        # ensure cftime
        var_ds = var_ds.convert_calendar(ds.time.dt.calendar, use_cftime=True)
        if allow_resample is True:
            var_ds = _adjust_frequency(var_ds, cf_freq, input_freq, time_cell_method)
        if rewrite_time_axis is True:
            var_ds = _rewrite_time_axis(var_ds, cf_freq)
        if "time" not in ds.cf.bounds and time_cell_method != "point":
            warn("adding time bounds")
            var_ds = _add_time_bounds(var_ds, cf_freq)
        if use_cftime is False:
            var_ds = var_ds.convert_calendar(ds.time.dt.calendar, use_cftime=False)
        var_ds = _set_time_encoding(var_ds, time_units, ds)

    if allow_units_convert is True:
        var_ds[out_name] = _cf_units_convert(
            var_ds[out_name], cmor_table, mapping_table
        )

    try:
        mapping = ds.cf["grid_mapping"]  # _get_pole(ds)
    except KeyError:
        warn(f"adding pole from archive specs: {domain_id}")
        mapping = _get_cordex_pole(domain_id)

    if "time" in mapping.coords:
        raise Exception("grid_mapping variable should have no time coordinate!")

    var_ds[mapping.name] = mapping

    return var_ds


def cmorize_variable(
    ds,
    out_name,
    cmor_table,
    dataset_table,
    mapping_table=None,
    grids_table=None,
    inpath=None,
    replace_coords=False,
    allow_units_convert=False,
    allow_resample=False,
    input_freq=None,
    domain_id=None,
    crop=None,
    time_units=None,
    rewrite_time_axis=False,
    outpath=None,
    **kwargs,
):
    """Cmorizes a variable.

    This functions call the python cmor API and creates a cmorized NetCDF file from the input
    dataset. Coordinates of the input dataset should be understandable by ``cf_xarray`` if they
    follow basic CF conventions. If a vertical coordinate is available, it should have an ``axis="Z"``
    attribute so it can be understood by ``cf_xarray`` and it should be named after the unique key of the
    coordinate in the cmor grids table (not the out_name), e.g., name it ``sdepth`` if you
    need a soil layer coordinate instead of ``depth``. Variables may have one of the following attributes
    that are used as keyword arguments in the call to `cmor_variable <https://cmor.llnl.gov/mydoc_cmor3_api/#cmor_variable>`_, e.g.,
    ``positive``, ``missing_value``, ``original_name``, ``history`` or ``comment``.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing at least the variable that should be cmorized.
    out_name: str
        CF out_name of the variable that should be cmorized. The corresponding variable name
        in the dataset is looked up from the mapping_table if provided.
    cmor_table : str or dict
        Cmor table dict of filepath to cmor table (json).
    dataset_table: str or dict
        Dataset table dict of filepath to dataset cmor table (json).
    mapping_table: dict
        Mapping of input variable names and meta data to CF out_name. Required if
        the variable name in the input dataset is not equal to out_name.
    grids_table: str
        Filepath to cmor grids table.
    inpath: str
        Path to cmor tables, if ``inpath == None``, inpath is the path
        to ``cmor_table``. This is required to find additional cmor tables,
        like ``CMIP6_coordinates``, ``CMIP6_grids`` etc.
    replace_coords: bool
        Replace coordinates from input file and create them from archive
        specifications. Only possible, if a domain identifier is given in the global
        attributes or as a keyword argument, e.g., see the ``domain_id`` keyword.
    allow_units_convert: bool
        Allow units to be converted if they do not agree with the
        units in the cmor table. Defaults to ``False`` to make the user aware of having
        correct ``units`` attributes set.
    allow_resample: bool
        Allow to resample temporal data to the frequency required from the cmor table.
        Handles both downsampling and upsampling. Defaults to ``False`` to make users aware
        of the correct frequency input.
    input_freq: str
        The frequency of the input dataset in pandas notation. It ``None`` and the dataset
        contains a time axis, the frequency will be determined automatically using
        ``pandas.infer_freq`` if possible.
    domain_id: str
        Cordex domain identifier. If ``None``, the domain will be determined by the ``domain_id``
        global attribute if available. If ``domain_id`` is given as a keyword, it will override
        a possible ``domain_id`` global attribute.
    crop: bool
        Crop dataset to official Cordex domain if it contains, e.g., a nudging zone. If set to ``None``,
        cropping will be done if a domain identifier is given or the domain can be identified automatically.
        Set to ``crop=False``, if you have an 'unofficial' domain_id.
    time_units: str
        Time units of the cmorized dataset (``ISO 8601``).
        If ``None``, time units will be set to default (``"days since 1950-01-01T00:00:00"``).
        If ``time_units='input'``, the original time units of the input dataset are used.
    rewrite_time_axis: bool
        Rewrite the time axis to CF compliant timestamps.
    outpath: str
        Root directory for output (can be either a relative or full path). This will override
        the outpath defined in the dataset cmor input table (``dataset_table``).
    **kwargs:
        Argumets passed to prepare_variable.

    Returns
    -------
    filename
        Filepath to cmorized file.

    Example
    -------

    To cmorize a dataset, you can use ,e.g.,::

        import cordex as cx
        from cordex.cmor import cmorize_variable
        from cordex.tables import cordex_cmor_table

        ds = cx.domain("EUR-11", dummy="topo").rename(topo="orog")
        dataset_table = cordex_cmor_table(f"CORDEX-CMIP6_remo_example")

        filename = cmorize_variable(
            ds,
            "orog",
            cmor_table=cordex_cmor_table("CORDEX-CMIP6_fx"),
            dataset_table=dataset_table,
            replace_coords=True,
            allow_units_convert=True,
        )

    """
    ds = ds.copy()

    if "CORDEX_domain" in kwargs:
        warn(
            "'CORDEX_domain' keyword is deprecated, please use the 'domain_id' keyword instead",
            DeprecationWarning,
            stacklevel=2,
        )
        domain_id = kwargs["CORDEX_domain"]
        del kwargs["CORDEX_domain"]

    if domain_id is None:
        try:
            domain_id = ds.cx.domain_id
        except Exception as e:
            warn(e)
            warn(
                "could not identify CORDEX domain, try to set the 'domain_id' if there is no domain identifier in global attributes."
            )
            crop = False
            replace_coords = False
    else:
        if "domain_id" in ds.attrs and ds.attrs.get("domain_id") != domain_id:
            warn(
                f"overwriting global attribute 'domain_id' with value '{ds.attrs['domain_id']}' with '{domain_id}' from keyword argument."
            )
        ds.attrs["domain_id"] = domain_id

    if domain_id and crop is None:
        crop = True

    if inpath is None:
        inpath = os.path.dirname(cmor_table)

    if not isinstance(dataset_table, dict):
        dataset_table = _read_table(dataset_table)

    if outpath:
        dataset_table["outpath"] = outpath

    ds_prep = prepare_variable(
        ds,
        out_name,
        cmor_table,
        domain_id=domain_id,
        mapping_table=mapping_table,
        replace_coords=replace_coords,
        input_freq=input_freq,
        rewrite_time_axis=rewrite_time_axis,
        time_units=time_units,
        allow_resample=allow_resample,
        allow_units_convert=allow_units_convert,
        crop=crop,
        **kwargs,
    )

    return cmorize_cmor(
        ds_prep, out_name, cmor_table, dataset_table, grids_table, inpath
    )


class CmorizerBase:
    def __init__(self):
        pass


class Cmorizer(CmorizerBase):
    def __init__(
        self,
        dataset_table=None,
        mapping_table=None,
        grids_table=None,
        inpath=None,
        replace_coords=False,
        allow_units_convert=False,
        allow_resample=False,
        input_freq=None,
        domain_id=None,
        time_units=None,
        rewrite_time_axis=False,
        outpath=None,
    ):
        """Cmorizes a variable.

        Parameters
        ----------
        dataset_table: str or dict
            Dataset table dict of filepath to dataset cmor table (json).
        mapping_table: dict
            Mapping of input variable names and meta data to CF out_name. Required if
            the variable name in the input dataset is not equal to out_name.
        inpath: str
            Path to cmor tables, if ``inpath == None``, inpath is the path
            to ``cmor_table``. This is required to find additional cmor tables,
            like ``CMIP6_coordinates``, ``CMIP6_grids`` etc.
        grids_table: str
            Filepath to cmor grids table.
        replace_coords: bool
            Replace coordinates from input file and create them from archive
            specifications.
        allow_units_convert: bool
            Allow units to be converted if they do not agree with the
            units in the cmor table. Defaults to ``False`` to make the user aware of having
            correct ``units`` attributes set.
        allow_resample: bool
            Allow to resample temporal data to the frequency required from the cmor table.
            Handles both downsampling and upsampling. Defaults to ``False`` to make users aware
            of the correct frequency input.
        input_freq: str
            The frequency of the input dataset in pandas notation. It ``None`` and the dataset
            contains a time axis, the frequency will be determined automatically using
            ``pandas.infer_freq`` if possible.
        domain_id: str
            Cordex domain short name. If ``None``, the domain will be determined by the ``domain_id``
            global attribute if available.
        time_units: str
            Time units of the cmorized dataset (``ISO 8601``).
            If ``None``, time units will be set to default (``"days since 1950-01-01T00:00:00"``).
            If ``time_units='input'``, the original time units of the input dataset are used.
        rewrite_time_axis: bool
            Rewrite the time axis to CF compliant timestamps.
        outpath: str
            Root directory for output (can be either a relative or full path). This will override
            the outpath defined in the dataset cmor input table (``dataset_table``).

        """
        super().__init__()
        self.dataset_table = dataset_table
        self.mapping_table = mapping_table
        self.grids_table = grids_table
        self.inpath = inpath
        self.replace_coords = replace_coords
        self.allow_units_convert = allow_units_convert
        self.allow_resample = allow_resample
        self.domain_id = domain_id
        self.time_units = time_units
        self.rewrite_time_axis = rewrite_time_axis
        self.outpath = outpath

        if op.isfile(self.dataset_table):
            self.dataset_table = _read_table(self.dataset_table)

        if outpath:
            self.dataset_table["outpath"] = outpath

    def preprocess(
        self,
        ds,
        out_name,
        cmor_table,
        replace_coords=False,
        allow_units_convert=False,
        allow_resample=False,
        input_freq=None,
        domain_id=None,
        time_units=None,
        rewrite_time_axis=False,
        use_cftime=False,
        squeeze=True,
    ):
        """prepares a variable for cmorization."""
        return prepare_variable(
            ds,
            out_name,
            cmor_table,
            mapping_table=self.mapping_table,
            replace_coords=replace_coords or self.replace_coords,
            allow_units_convert=allow_units_convert or self.allow_units_convert,
            allow_resample=allow_resample or self.allow_resample,
            input_freq=input_freq,
            domain_id=domain_id or self.domain_id,
            time_units=time_units or self.time_units,
            rewrite_time_axis=rewrite_time_axis or self.rewrite_time_axis,
            use_cftime=use_cftime,
            squeeze=squeeze,
        )

    def _write_with_cmor(self, ds, out_name, cmor_table):
        return cmorize_cmor(
            ds, out_name, cmor_table, self.dataset_table, self.grids_table, self.inpath
        )

    def cmorize(self, ds, out_name, cmor_table):
        """Cmorizes a variable.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset containing at least the variable that should be cmorized.
        out_name: str
            CF out_name of the variable that should be cmorized. The corresponding variable name
            in the dataset is looked up from the mapping_table if provided.
        cmor_table : str or dict
            Cmor table dict of filepath to cmor table (json).
        dataset_table: str or dict
            Dataset table dict of filepath to dataset cmor table (json).
        mapping_table: dict
            Mapping of input variable names and meta data to CF out_name. Required if
            the variable name in the input dataset is not equal to out_name.
        grids_table: str
            Filepath to cmor grids table.
        inpath: str
            Path to cmor tables, if ``inpath == None``, inpath is the path
            to ``cmor_table``. This is required to find additional cmor tables,
            like ``CMIP6_coordinates``, ``CMIP6_grids`` etc.
        replace_coords: bool
            Replace coordinates from input file and create them from archive
            specifications.
        allow_units_convert: bool
            Allow units to be converted if they do not agree with the
            units in the cmor table. Defaults to ``False`` to make the user aware of having
            correct ``units`` attributes set.
        allow_resample: bool
            Allow to resample temporal data to the frequency required from the cmor table.
            Handles both downsampling and upsampling. Defaults to ``False`` to make users aware
            of the correct frequency input.
        input_freq: str
            The frequency of the input dataset in pandas notation. It ``None`` and the dataset
            contains a time axis, the frequency will be determined automatically using
            ``pandas.infer_freq`` if possible.
        domain_id: str
            Cordex domain short name. If ``None``, the domain will be determined by the ``domain_id``
            global attribute if available.
        time_units: str
            Time units of the cmorized dataset (``ISO 8601``).
            If ``None``, time units will be set to default (``"days since 1950-01-01T00:00:00"``).
            If ``time_units='input'``, the original time units of the input dataset are used.
        rewrite_time_axis: bool
            Rewrite the time axis to CF compliant timestamps.
        outpath: str
            Root directory for output (can be either a relative or full path). This will override
            the outpath defined in the dataset cmor input table (``dataset_table``).
        **kwargs:
            Argumets passed to prepare_variable.

        Returns
        -------
        filename
            Filepath to cmorized file.

        """
        ds_prep = self.preprocess(ds, out_name, cmor_table)
        return self._write_with_cmor(ds_prep, out_name, cmor_table)

    def cfinfo(self, out_name, cmor_table):
        if isinstance(cmor_table, str):
            cmor_table = _read_table(cmor_table)
        return _get_cfvarinfo(out_name, cmor_table)
