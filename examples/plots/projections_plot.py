"""projections plot example script
"""

import xarray as xr
from cordex.plot import projections as cxproj

infile ='data/remo_emhh_oro.nc'
varname='var129'

#Example for proj_Plot(proj) using different projections 
import cartopy.crs as ccrs  
projections = [ccrs.RotatedPole(),ccrs.PlateCarree(),
               ccrs.Robinson(),
               ccrs.Mercator(),
               ccrs.Orthographic(),ccrs.NearsidePerspective()]

for proj in projections:
    cxproj.proj_Plot(infile, varname, proj)

#Global plots
for proj in projections:
    cxproj.proj_Plot(infile, varname, proj,'global')

#Specifying projections coordinates
cxproj.proj_Plot(infile, varname,ccrs.Orthographic(10,45),'(10,45)_global')