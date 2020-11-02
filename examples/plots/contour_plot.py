"""contour plot example script
"""

import xarray as xr
import cordex.plot.contour as cxplt

infile ='data/remo_emhh_oro.nc'
#ds1 = xr.open_dataset(infile1)
#ds1.info()

# call plotting routine
varname='var129'

cxplt.contour2(infile, varname, 'Orographie_EU')
cxplt.contour3(infile, varname, 'Orographie_EU_Subplots')


# The other domains: Europe, Africa, North and South America, Asia
my_oro_rot = "data/remo_europe_r_oro.nc"
af_oro_rot = "data/remo_africa_r_oro.nc"    
na_oro_rot = "data/remo_n_america_r_oro.nc"
sa_oro_rot = "data/remo_s_america_r_oro.nc"
as_oro_rot = "data/remo_s_asia_r_oro.nc"
infile1=[my_oro_rot,af_oro_rot,na_oro_rot,sa_oro_rot,as_oro_rot]
cxplt.multi_domain(infile1, varname, 'World')

# isotherm contour lines over Europe (two subplots with contour lines and filled contours)
infile2='data/remo_emhh_temp.nc'
varname2='var167' 
cxplt.contour7(infile2, varname2, 'Temp. filled Contours in [ºC]')

# Orography with isotherm contour lines over Europe
infile1 ='data/remo_emhh_oro.nc'
infile2 = "data/remo_emhh_temp.nc"
varname1='var129'
varname2='var167'
cxplt.contour8(infile1,infile2,varname1,varname2, 'Surface Geopotential in [m] and Temperature in [K]')

