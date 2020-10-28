"""labeling, coloring, leveling subroutines
"""

import numpy as np
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def level_color_label(varname):
    """*Function*
    Sets the colormapping, leveling and labeling of the variable chosen.
    **Arguments:**
        *varname:*
            Name of the variable to read.

    **Returns:**
        *clevs*
            Levels or values at which the variable is depicted.
        *col*
            Colormap with which the variable is shown.
        *label*
            Gives the units and name of the variable.
    """
    
    #Orographie
    if varname == 'var129':
        # define custom levels and colormap in [m]
        clevs = [[1, 25, 50, 100, 150, 200, 300, 500, 750,
             1000, 1250, 1500, 1750, 2000, 2500, 3000],[1, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000,
             3500, 4000, 4500, 5000]]
        colmap = np.array([[0.39, 0.58, 0.93], [0.00, 0.39, 0.00],[0.16, 0.49, 0.00], [0.31, 0.59, 0.00],
                       [0.47, 0.69, 0.00], [0.63, 0.78, 0.00],[0.78, 1.00, 0.00], [1.00, 1.00, 0.00],
                       [0.90, 0.86, 0.00], [0.78, 0.71, 0.00],[0.67, 0.59, 0.00], [0.57, 0.43, 0.00],
                       [0.47, 0.29, 0.00], [0.35, 0.16, 0.00],[0.53, 0.43, 0.35], [0.71, 0.71, 0.71],
                       [1.00, 0.98, 0.98]], 'f')
        col=ListedColormap(colmap)
        label='Surface Geopotential [m]'
    #Temperature in [ºC]
    if varname == 'var167':
        clevs=[-30., -27., -24., -21., -18., -15., -12., -9., -6., -3., 0., 3., 6., 9., 12., 15., 18.,22.]
        col='rainbow'
        label='Temperature [ºC]'
    #Pressure in [hPa]
    if varname == 'var151':
        clevs = [1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015,
             1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025]
        col='jet'
        label='Pressure [hPa]'
    #Streamvector strength
    if varname == 'var171':
        clevs = np.arange(0, 10, 0.5)
        col='jet'
        label='Vector Length [m]'
    if varname == 'var165':
        clevs = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        col = 'jet'
        label = '10m u-velocity [$\frac{m}{s}$]'
    if varname == 'var166':
        clevs = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        col = 'jet'
        label = '10m u-velocity [$\frac{m}{s}$]'  
    
    return clevs, col, label