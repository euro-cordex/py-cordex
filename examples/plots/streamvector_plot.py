"""contour plot example script
"""

import xarray as xr
import cordex.plot.stream_vector as cxsv

infile="data/remo_wind.nc"
var1='var165' 
var2='var166'
var3='var151'
var4='var171'

#cxsv.stream_and_vector1(infile, var1,var2,var3,"Stream and Pressure") 

#cxsv.stream_and_vector2(infile, var1,var2,var3, "Vector and Pressure") 

cxsv.stream_and_vector3(infile, var1,var2,var3,var4,"10m wind speed")
 
cxsv.stream_and_vector4(infile, var1,var2,"Wind speed with vector length meassure")

cxsv.stream_and_vector5(infile, var1,var2,var3,"Mean Sea Level Pressure in [Pa]") 



