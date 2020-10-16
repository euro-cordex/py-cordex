"""plot example script
"""

import xarray as xr
import cordex.plot as cxplt


infile1 ='data/remo_emhh_oro.nc'
ds1 = xr.open_dataset(infile1)
ds1.info()

# call plotting routine
cxplt.contour2(infile1, 'var129', 'Orographie_EU')
