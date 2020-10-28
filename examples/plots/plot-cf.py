"""plot example script
"""

import xarray as xr
from netCDF4 import Dataset
import cordex.plot as cxplt


from cordex.plot import grid as tx


infile1 ='data/orog_EUR-11_ECMWF-ERAINT_evaluation_r0i0p0_GERICS-REMO2015_v1_fx.nc'
infile2 ='data/remo_emhh_oro.nc'

#ds1 = xr.open_dataset(infile1)

infile = infile2

#ds1.info()

# call plotting routine

orog_xr = xr.open_dataset(infile) #['orog']
orog_ds = Dataset(infile)

mapping = tx.get_grid_mapping(orog_xr, 'var129')
print(mapping)
print(mapping.grid_north_pole_latitude)


mapping = tx.get_grid_mapping(orog_ds, 'var129')
print(mapping)
print(mapping.grid_north_pole_latitude)

cxplt.contour2(infile, 'var129', 'Orographie_EU')

#print(mapping)


