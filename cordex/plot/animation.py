"""plotting subroutines
"""

import cordex.plot.read_nc as read_nc
import cordex.plot.level_color_label as level_color_label
import matplotlib.pyplot as plt
import numpy as np


try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,LatitudeLocator)
    import cartopy.feature as cfeature
except:
    print('cartopy not installed, plotting capabilities reduced...')



def level_color_label_ERA5(var,steps,*args):
    """*Function*
    Sets the colormapping, leveling and labeling of the variable chosen. If there are
    more variables to be ploted, fill in the optional arguments (see example Plot5 with multi_domain())
    **Arguments:**
        *varname:*
            Variable to be plotted. Automatic input by function read_nc().
        *steps:*
            Steps for the colorbar.
        *args:*
            Optional argument, additional varibale to be plotted.
    **Returns:**
        *clevs*
            Levels or values at which the variable is depicted. Array from min to max value with integers.
        *col*
            Colormap with which the variable is shown.
        *label*
            Gives the units and name of the variable.
    """
    var1=np.append(var,args)
    clevs=np.linspace(np.floor(var1.min()), np.ceil(var1.max()), steps,dtype=int)
    col='rainbow'
    label=var.long_name +' ['+ var.units+']'
    return clevs, col, label  

def contour_video(infile, varname):
    """*Video Start 1*
    Gives a snapshot of the evolution of the colormap as contour plot of the given 
    file of the given variable to later on make a video.
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.
    **Returns:**
        *plot*
            Returns a snapshot of the chosen variable and file to produce the video.
    """
    #read netcdf and select variable
    var,rot_pole = read_nc(infile, varname)
    
    #levels, colors, label 
    clevs,col,label= level_color_label_ERA5(var,20)
    
    ax = plt.axes(projection=ccrs.PlateCarree())

    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #Country boarders
    ax.add_feature(cfeature.BORDERS)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)

    #Actual Plot
    cax=var[0,:,:].plot(ax=ax,levels=clevs, transform=ccrs.PlateCarree(),animated=True,cmap=col,
             cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=1, label=label),infer_intervals=True)
    
    return cax,var
