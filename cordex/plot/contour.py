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

def contour2(infile, varname,title):
    """*Plot 2*
    Contour-plot with user defined colors, contour levels and as
    pdf output. 
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.
        *title:*
            Title of the figure as string.
    **Returns:**
        *plot*
            Returns a pdf with the contour-plot.
    """
    #read netcdf and select variable
    var,rot_pole = read_nc(infile, varname)
    #levels, colors, label 
    clevs,col,label= level_color_label(varname)
    
    fig = plt.figure(figsize=(10,8))
    
    #Rotated pole projection with Cartopy
    rotated_pole = ccrs.RotatedPole(pole_latitude=rot_pole[0], 
                                pole_longitude=rot_pole[1])
    ax = plt.axes(projection=rotated_pole)
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #Country boarders
    ax.add_feature(cfeature.BORDERS)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)

    #Actual Plot
    var.plot(ax=ax,vmin=0,transform=rotated_pole,cmap=col,levels=clevs[0],
             cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=1, label=label))
    ax.set_title(title)
    
    fig.savefig('results/'+title+'.pdf')


def contour3(infile, varname, title):
    """*Plot 3*
    Two subplots as contour-plots with user defined colors, contour levels and as
    pdf output.
    Subplot 1: Map of the Globus showing the region of interest.
    Subplot2: Only region of interest
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.
        *title:*
            Title of the figure as string.
    **Returns:**
        *plot*
            Returns a pdf with the contour-plot.
    """
    #read netcdf and select variable
    var,rot_pole = read_nc(infile, varname)
    #levels, colors, label 
    clevs,col,label= level_color_label(varname)
    
    #Setting the subplots
    fig = plt.figure(figsize=(13,10))

    #rotated pole
    rotated_pole = ccrs.RotatedPole(pole_latitude=rot_pole[0], 
                                pole_longitude=rot_pole[1])
    #Rotated pole transform and satellite projection with Cartopy with default defined colormap
    ax1 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(10,45))
    ax1.coastlines(resolution='50m', color='black', linewidth=1)
    ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    var.plot(ax=ax1,vmin=0,transform=rotated_pole,cmap=col,levels=clevs[0],
             cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=0.5, label=label))
    #Country boarders
    ax1.add_feature(cfeature.BORDERS)
    ax1.set_title(title+" from Satellite")
    ax1.set_global()

    #Rotated pole projection with Cartopy, with half of the default defined colormap
    ax2 = plt.subplot(1, 2, 2, projection=rotated_pole)
    #printing coastlines
    ax2.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    #Actual Plot
    var.plot(ax=ax2,vmin=0,transform=rotated_pole,cmap=col,levels=clevs[0][::2],
             cbar_kwargs=dict(orientation='horizontal',pad=0.1, shrink=1, label=label))
    #Country boarders
    ax2.add_feature(cfeature.BORDERS)
    ax2.set_title(title)
    fig.savefig('results/'+title+'.pdf')

def multi_domain(infile, varname,title):
    """*Plot 5*
    Ploting different regions of interest in one map of the globus with rotated pole coordinates for each region. 
    !! infile is a list/ array of datasets
    **Arguments:**
        *infile:*
            Listed path of each File/ Netcdf to read. Must be an array or list.
        *varname:*
            Name of the variable to read.
        *title:*
            Title of the plot.
    **Returns:**
        *plot*
            Returns the contour-plot of the globus with all regions of interest.
    """
    fig = plt.figure(figsize=(15,12))
    
    for i in infile:
    #read netcdf and select variable
        var,rot_pole = read_nc(i, varname)
        #levels, colors, label 
        clevs,col,label= level_color_label(varname)
        rotated_pole = ccrs.RotatedPole(pole_latitude=rot_pole[0], 
                                pole_longitude=rot_pole[1])
        ax = plt.axes(projection=ccrs.Robinson())

    #Actual Plot
        var.plot(vmin=0,transform=rotated_pole,cmap=col,levels=clevs[1],add_colorbar=False)
        ax.set_global()
    
    var.plot(ax=ax,vmin=0,transform=rotated_pole,cmap=col,levels=clevs[1],
             cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=0.5, label=label))
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #Country boarders
    ax.add_feature(cfeature.BORDERS)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
       
    ax.set_title(title)
    fig.savefig('results/'+title+'.pdf')



def contour7(infile, varname, title):
    """*Plot 7*
    Two subplots as contour-plots 
    Plot1 shows only temp. contour lines with values in ºC with rainbow colomapping. 
    Plot2 filled contours with colormapping 'YlOrBr'.
    **Arguments:**
        *infile:*
            Path of the File/ Netcdf to read.
        *varname:*
            Name of the variable to read.

    **Returns:**
        *plot*
            Returns a pdf with the contour-plot.
    """
    #read netcdf and select variable
    var,rot_pole = read_nc(infile, varname)
    #levels, colors, label 
    clevs,col,label= level_color_label(varname)
    var=var-275.
    #Setting the subplots
    fig = plt.figure(figsize=(13,10))
    clevs1=[-30., -27., -24., -21., -18., -15., -12., -9., -6., -3., 0., 3., 6., 9., 12., 15., 18.,22.]


    #rotated pole
    rotated_pole = ccrs.RotatedPole(pole_latitude=rot_pole[0], pole_longitude=rot_pole[1])
    #Rotated pole transform and satellite projection with Cartopy with default defined colormap
    ax1 = plt.subplot(1, 2, 1, projection=rotated_pole)
    ax1.coastlines(resolution='50m', color='black', linewidth=1)
    ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)    
    a1=ax1.contour(var.rlon, var.rlat, var, levels=clevs[::2],cmap='rainbow',
             transform=rotated_pole)
    ax1.clabel(a1, a1.levels, fmt='%1.f', inline=True, fontsize=10,use_clabeltext=True)

    #Country boarders
    ax1.add_feature(cfeature.BORDERS)
    ax1.set_title("Temp. Contours in [ºC]")


    #Rotated pole projection with Cartopy, with half of the default defined colormap
    ax2 = plt.subplot(1, 2, 2, projection=rotated_pole)
    #printing coastlines
    ax2.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)
    #Actual Plot
    a2=ax2.contourf(var.rlon, var.rlat, var,levels=clevs,cmap='YlOrBr',transform=rotated_pole)
    fig.colorbar(a2, ax=ax2,shrink=0.5)
    #Country boarders
    ax2.add_feature(cfeature.BORDERS)
    ax2.set_title(title)
    fig.savefig('results/'+title+'.pdf')

def contour8(infile1, infile2, varname1, varname2,title):
    """*Plot 8*
    The output is an image with the map of Europe with the orographie in the
    background and temperature as contour lines entitled
    *2m temperature* because we overlayed 2m temperature as the last variable
    before plotting.
    **Arguments:**
        *infile1:*
            Path of the first File/ Netcdf to read, in this case orography
        *infile2:*
            Path of the second File/ Netcdf to read, in this case temperature
        *varname1:*
            Name of the first variable to read.
        *varname2:*
            Name of the second variable to read.
    **Returns:**
        *plot*
    """  
    #read netcdf and select variable
    var1,rot_pole1 = read_nc(infile1, varname1)
    #levels, colors, label 
    clevs1,col1,label1= level_color_label(varname1)
    
    fig = plt.figure(figsize=(10,8))
    
    #Rotated pole projection with Cartopy
    rotated_pole1 = ccrs.RotatedPole(pole_latitude=rot_pole1[0], 
                                pole_longitude=rot_pole1[1])
    ax = plt.axes(projection=rotated_pole1)
    #printing coastlines
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #setting gridlines
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, y_inline=False,x_inline=False)

    #Actual Plot
    var1.plot(ax=ax,vmin=0,transform=rotated_pole1,cmap=col1,levels=clevs1[0],
             cbar_kwargs=dict(orientation='vertical',pad=0.1, shrink=1, label=label1))
    
    ##### Temperature Contours
    #read netcdf and select variable
    var2,rot_pole2 = read_nc(infile2, varname2)
    var2=var2-275.
    #levels, colors, label 
    clevs2,col2,label2= level_color_label(varname2)
    a=ax.contour(var2.rlon, var2.rlat, var2, levels=clevs2[::2],colors='k',
             transform=rotated_pole1)
    clabels=ax.clabel(a, a.levels, fmt='%1.f', inline=True, fontsize=10,use_clabeltext=True)
    ax.set_title(title)
    #draw country boarders
    ax.add_feature(cfeature.BORDERS)

    fig.savefig('results/'+title+'.pdf')
