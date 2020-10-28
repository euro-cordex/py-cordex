"""plotting subroutines
"""

import cordex.plot.read_nc as read_nc
import cordex.plot.level_color_label as level_color_label
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np

try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,LatitudeLocator)
    import cartopy.feature as cfeature
except:
    print('cartopy not installed, plotting capabilities reduced...')



def stream_and_vector1(infile, varname1,varname2,varname3,title):
    """*Streamplot 1*
        Mean-sea-level pressure with wind streamlines as overlay on a
        non-rotated grid. 
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname1:*
            Name of the first variable to read, stream in x-direction
        *varname2:*
            Name of the second variable to read, stream in y-direction
        *varname3:*
            Name of the third variable to read, pressure
    **Returns:**
        *plot*
        Stream plot with pressure colormap
    """  
    #read netcdf and select variable
    var1,rot_pole1 = read_nc(infile, varname1)
    var2,rot_pole2 = read_nc(infile, varname2)
    var3,rot_pole3 = read_nc(infile, varname3)
    
    fig = plt.figure(figsize=(10,8))
    
    # convert the pressure to hPa
    var3 = var3 / 100.
    
    clevs = [1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
             1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025]
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    
    a=var3.plot.contourf(ax=ax, levels=clevs, cmap = "jet", transform=ccrs.PlateCarree())

    ax.streamplot(var1.lon,var2.lat,var1.data,var2.data, 
                  linewidth=1, density=1, color='k')

    ax.set_title(title)
    #draw country boarders
    ax.add_feature(cfeature.BORDERS)
    fig.savefig('results/'+title+'.pdf')


def stream_and_vector2(infile, varname1,varname2,varname3, title):
    """*Streamplot 2*
         Mean-sea-level pressure with wind vectors as overlay on a non-rotated grid. 
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname1:*
            Name of the first variable to read, stream in x-direction
        *varname2:*
            Name of the second variable to read, stream in y-direction
        *varname3:*
            Name of the third variable to read, pressure
    **Returns:**
        *plot*
        Vector plot with pressure colormap
    """  
    #read netcdf and select variable
    var1,rot_pole1 = read_nc(infile, varname1)
    var2,rot_pole2 = read_nc(infile, varname2)
    var3,rot_pole3 = read_nc(infile, varname3)
    
    fig = plt.figure(figsize=(10,8))
    
    # convert the pressure to hPa
    var3 = var3 / 100.
    
    clevs = [1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
             1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025]
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    
    a=var3.plot.contourf(ax=ax, levels=clevs, cmap = "jet",shrink=0.5, transform=ccrs.PlateCarree(),add_colorbar=False)
    cbar = plt.colorbar(a,ax=ax,shrink=0.9,label='var151')
    
    Q = plt.quiver( var1.lon,var2.lat,var1.data,var2.data, 
               pivot='middle', 
               transform=ccrs.PlateCarree(), 
               regrid_shape=15,
               angles='xy', scale_units='xy', scale=1)

    # plot quiver key
    qk = plt.quiverkey(Q, 
                   1.07, 0.99,                  # x,y label position
                   1., r'$1 \frac{m}{s}$', # choose units + update string
                   labelpos='E',                # add label to the right
                   coordinates='axes'
                   )   
    
    ax.set_title(title)
    #draw country boarders
    ax.add_feature(cfeature.BORDERS)
    fig.savefig('results/'+title+'.pdf')


def stream_and_vector3(infile, varname1,varname2,varname3,varname4, title):
    """*Streamplot 3*
         Wind direction as streamlines with windspeed as color and costum streamline
    level spacing.
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname1:*
            Name of the first variable to read, stream in x-direction
        *varname2:*
            Name of the second variable to read, stream in y-direction
        *varname3:*
            Name of the third variable to read, pressure
        *varname4:*
            Name of the third variable to read, strength
    **Returns:**
        *plot*
        Wind direction as streamlines with windspeed as color and costum streamline
    level spacing.
    """  
    #read netcdf and select variable
    var1,rot_pole = read_nc(infile, varname1)
    var2,rot_pole = read_nc(infile, varname2)
    var3,rot_pole = read_nc(infile, varname3)
    var4,rot_pole = read_nc(infile, varname4)

    fig = plt.figure(figsize=(10,8))
    
    # convert the pressure to hPa
    var3 = var3 / 100.
    
    clevs = [1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
             1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025]
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    
    #lenght of the velocity vectors
    length= (var1.data**2 +var2.data**2)**0.5
    stream=ax.streamplot(var1.lon,var2.lat,var1.values,var2.values, 
                  linewidth=1, density=30,color=var4.values, cmap='jet')
    colorbounds = np.arange(0, 11, 1)
    cb = fig.colorbar(stream.lines,ax=ax,shrink=0.9,label='10m wind speed in [m/s]',
                  boundaries=colorbounds, 
                  ticks=colorbounds,
                  spacing='uniform',
                  orientation='vertical')
    
    ax.set_title(title)
    fig.savefig('results/'+title+'.pdf')


def stream_and_vector4(infile, varname1,varname2, title):
    """*Streamplot 4*
         Wind direction as streamlines with windspeed as color and costum streamline
    level spacing.
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname1:*
            Name of the first variable to read, stream in x-direction
        *varname2:*
            Name of the second variable to read, stream in y-direction
    **Returns:**
        *plot*
        Wind direction as streamlines with windspeed as color and costum streamline
    level spacing.
    """  
    #read netcdf and select variable
    var1,rot_pole = read_nc(infile, varname1)
    var2,rot_pole = read_nc(infile, varname2)
    fig = plt.figure(figsize=(10,8))

    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    
    #vector length
    length= (var1.data**2 +var2.data**2)**0.5

    Q = ax.quiver( var1.lon,var2.lat,var1.data,var2.data,length,
               pivot='middle', 
               transform=ccrs.PlateCarree(), 
               cmap = "jet" ,angles='xy', scale_units='xy', scale=1,regrid_shape=100  )
    ax.set_extent([-10, 30, -40, 10], ccrs.PlateCarree())
    # plot quiver key
    qk = ax.quiverkey(Q, 
                   1.07, 0.99,                  # x,y label position
                   1., r'$1 \frac{m}{s}$', # choose units + update string
                   labelpos='E',                # add label to the right
                   coordinates='axes'
                   )
    colorbounds = np.arange(0, 10, 0.5)
    cbar = plt.colorbar(Q,
                    ax=ax,
                    boundaries=colorbounds, 
                    ticks=colorbounds,
                    orientation='vertical',shrink=0.7, label='Vector length [m]')
    ax.set_title(title)
    #draw country boarders
    ax.add_feature(cfeature.BORDERS)
    fig.savefig('results/'+title+'.pdf')


def stream_and_vector5(infile, varname1,varname2,varname3,title):
    """*Streamplot 5*
         Only the vectors from plot 2 with coloring of the vectors depending
    on the vector pressure. 
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname1:*
            Name of the first variable to read, stream in x-direction
        *varname2:*
            Name of the second variable to read, stream in y-direction
        *varname3:*
            Name of the third variable to read, pressure
    **Returns:**
        *plot*
        Only the vectors from plot 2 with coloring of the vectors depending
    on the vector pressure. 
    """  
    #read netcdf and select variable
    var1,rot_pole = read_nc(infile, varname1)
    var2,rot_pole = read_nc(infile, varname2)
    var3,rot_pole = read_nc(infile, varname3)
    fig = plt.figure(figsize=(10,8))

    # convert the pressure to hPa
    var3 = var3 / 100.
    
    clevs = [1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
             1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025]
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    
    Q = ax.quiver( var1.lon,var2.lat,var1.data,var2.data,var3.data,
               pivot='middle', 
               transform=ccrs.PlateCarree(), 
               cmap = "jet" ,angles='xy', scale_units='xy', scale=1,regrid_shape=100)
    ax.set_extent([-10, 28, -40, 10], ccrs.PlateCarree())
    # plot quiver key
    qk = ax.quiverkey(Q, 
                   1.07, 0.99,                  # x,y label position
                   1., r'$1 \frac{m}{s}$', # choose units + update string
                   labelpos='E',                # add label to the right
                   coordinates='axes'
                   )
    colorbounds = np.arange(0, 10, 0.5)
    cbar = plt.colorbar(Q,
                    ax=ax,
                    boundaries=clevs, 
                    ticks=clevs,
                    orientation='vertical',shrink=0.7, label='Pressure [Pa]')
    ax.set_title(title)
    #draw country boarders
    ax.add_feature(cfeature.BORDERS)
    fig.savefig('results/'+title+'.pdf')
