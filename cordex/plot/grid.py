

import xarray as xr



GRID_MAPPING = 'grid_mapping'


def get_grid_mapping(ds, varname=None):
    """Returns the grid mapping variable of a netcdf variable or dataset.
    """
    if varname:
        return get_grid_mapping_from_variable(ds, varname)
    else:
        return get_grid_mapping_from_dataset(ds)
#    if hasattr(da, GRID_MAPPING):


def get_grid_mapping_from_variable(ds, varname):
    """Returns the grid mapping variable associated with a netcdf variable.
    """
    mapping_name = get_grid_mapping_name_from_variable(ds[varname])
    if mapping_name is not None and mapping_name in ds.variables:
        return ds[mapping_name]
    else:
        return None

def get_grid_mapping_from_dataset(ds):
    """Returns the grid mapping variables from a dataset.
    """
    return [mapping for mapping in (get_grid_mapping(ds, var) for var in ds.variables) if mapping is not None]


def get_grid_mapping_name_from_variable(var):
    try:
        return getattr(var, GRID_MAPPING)
    except:
        return None


#def get_grid_mapping_names_from_dataset(ds):
#    return [mapping for mapping in (get_grid_mapping_name_from_variable(ds[var]) for var in ds.variables) if mapping is not None]
