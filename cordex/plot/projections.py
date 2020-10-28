"""plotting subroutines
"""

import cordex.plot.read_nc as read_nc
import cordex.plot.level_color_label as level_color_label
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,LatitudeLocator)
    import cartopy.feature as cfeature
except:
    print('cartopy not installed, plotting capabilities reduced...')

def proj_Plot(infile, varname, proj, *args):
    """*Plot 4*
    Ploting different projections. For the map fo the globus, uncomment "ax.set_global()".
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.
        *proj:*
            Wanted cartopy projection (e.g ccrs.RotatedPole(),ccrs.PlateCarree(),ccrs.Robinson(),
               ccrs.Mercator(),ccrs.Orthographic(),ccrs.NearsidePerspective())
        *args:* If there is an argument, the plot is of the whole Globus, for ex. 'global'.
                No entry, delivers plots of the specific region.
    **Returns:**
        *plot*
            Returns the contour-plot with the desired projection.
    """
    fig = plt.figure(figsize=(10,8))
    
    #read netcdf and select variable
    var,rot_pole = read_nc(infile, varname)
    #levels, colors, label 
    clevs,col,label= level_color_label(varname)
    
    if proj==ccrs.RotatedPole():
        proj = ccrs.RotatedPole(pole_latitude=rot_pole[0], 
                                pole_longitude=rot_pole[1])        
    else:
        pass
    rotated_pole = ccrs.RotatedPole(pole_latitude=rot_pole[0], pole_longitude=rot_pole[1])  

    ax = plt.axes(projection=proj)
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    import cartopy.feature as cfeature
    ax.add_feature(cfeature.BORDERS)
    var.plot(ax=ax,vmin=0,transform=rotated_pole,cmap=col,levels=clevs[0],cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=0.5, label=label))
    ax.set_title(f'{type(proj)}')
    if args:
        ax.set_global()
    else:
        pass

    fig.savefig('results/'+f'{type(proj)}'[20:-2]+'_'+str(args)[2:-3]+'.pdf')

    
