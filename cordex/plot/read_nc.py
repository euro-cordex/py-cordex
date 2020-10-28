"""netcdf reading subroutines
"""

import xarray as xr

#Read NetCDF
def read_nc(infile, varname):
    """Read NetCDF
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.

    **Returns:**
        *var:*
            Variable to be read
        *rot_pole:*
            Coordinates of rotated pole if there is a rotated pole defined in the 
            Netcdf in form [rot_latitude,rot_longitude].
    """
    ds1 = xr.open_dataset(infile)
    var = ds1.__getitem__(varname).squeeze()
    try: 
        rot_pole = [rot_latitude,rot_longitude] = [ds1.rotated_pole.grid_north_pole_latitude, ds1.rotated_pole.grid_north_pole_longitude]
        return var, rot_pole
    except:
        pass
        return var, None
#Example Africa without rotated pole -> Returns only var
#var,rot_pole = read_nc(my_wind,'var165')
#Example Europe with rotated pole -> Returns var and rot_pole
#var,rot_pole = read_nc(infile1,'var129')