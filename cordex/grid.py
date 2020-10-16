# -*- coding: utf-8 -*-
# flake8: noqa
"""
Defining a grid for rotated and non-rotated coordinates.
RotGrid is outdated because of Code duplication.
"""
import logging
import numpy as np
import math

from cordex import __version__

__author__ = "Kevin Sieck, Claas Teichmann, Lars Buntemeyer"
__copyright__ = "Lars Buntemeyer"
__license__ = "mit"

_logger = logging.getLogger(__name__)


class Grid(object):
    """This class contains gridded geographic coordinates.

    Variables ending with _geo are describing the points in the geographical
    non-rotated coordinates system you find on a globe.

    **Attributes:**
        *lon_arr:*
            longitudes in geographical system (2d-array)
        *lat_arr:*
            latitudes in geographical system (2d-array)
        *cyclic:*
            True if the data is on a cyclic (global) grid


    Last changes 12.08.2019 by Lars Buntemeyer

    .. note::
        Grid class is now generalized for holding arbitrary grid info,
        a 'rotated grid' is now simply one with pol_lat != 90. Default
        is not rotated (pol_lat = 90.). Ann.: a grid instance should only
        hold one valid dataset. This dataset can be transformed by creating
        a new grid instance, e.g., my_grid.transform().

    """


    # Methods
    def __init__(self, lon_arr, lat_arr, pol_lon=None, pol_lat=None):
        """Setting lon/lat-array

        **Arguments:**
            *lon_arr:*
                longitudes in geographical coordinate system (2d-array)
            *lat_arr:*
                latitudes in geographical coordinate system (2d-array)
            *pol_lon:*
                longitude of North Pole (Default: 180, not rotated)
            *pol_lat:*
                latitude of North Pole (Default: 90, not rotated)
        """
        self.pol_lon = 180. if pol_lon is None else pol_lon
        self.pol_lat =  90. if pol_lat is None else pol_lat
        self.lon_arr_geo, self.lat_arr_geo = self.init_lon_lat_arr(lon_arr, lat_arr)
        self.lon_arr, self.lat_arr         = self.init_lon_lat_arr(lon_arr, lat_arr)
        assert(self.lon_arr_geo.shape == self.lat_arr_geo.shape)
        assert(self.lon_arr.shape == self.lat_arr.shape)
        self._check_cyclic(self.lon_arr_geo)


    def __eq__(self, other):
        """Check for equality.

        Two objects are equal if the grids are equal.
        """
        if self.get_dimensions() == other.get_dimensions():
            is_equal = (np.allclose(self.lon_arr, other.lon_arr) and
                        np.allclose(self.lat_arr, other.lat_arr))
        else:
            is_equal = False
        return is_equal


    def __str__(self):
        """Details of the Grid to a string.
        """
        out_tmplt = (
                "Pole (lon/lat): {pollon}/{pollat}\n"
                "lon_arr:\n{lon_arr}\n"
                "lat_arr:\n{lat_arr}\n"
                )
        dic = {'pollon': self.pol_lon,
               'pollat': self.pol_lat,
               'lon_arr': self.lon_arr,
               'lat_arr': self.lat_arr
               }
        return out_tmplt.format(**dic)


    @property
    def rotated(self):
        """Checks if the pole is rotated.

        We ignore the longitude and simply check if the latitude is 90 degrees.
        """
        return self.pol_lat != 90.


    def _check_cyclic(self, lon_arr):
        """Sets the *cyclic* attribute.

        Checks if a grid is cylic, which means that the longitude span
        is 360 degree. If it is cyclic the *cyclic* attribute will be set
        to *True*.

        **Arguments:**
            *lon_arr:*
                The longitude array of the grid.

        """
        lon_step = abs(lon_arr[0, 0] - lon_arr[0, 1])
        xdim = lon_arr.shape[1]
        if xdim * lon_step < 360:
            self.cyclic = False
        else:
            self.cyclic = True


    def is_cyclic(self):
        """Returns *True* if the grid is cyclic.

        **Returns:**
            *cyclic:*
                *True* if the grid is cyclic. Otherwise *False*.
        """
        return self.cyclic


    @staticmethod
    def init_lon_lat_arr(lon_arr, lat_arr):
        """Creates two 2d-arrays of two 1d-lon/lat-arrays

        Two 2d-arrays are created if the dimension of both input-arrays is one,
        otherwise the arrays are just returned.
        Values are filled up to the size of the other array to create 2d-arrays
        with dimensions:
        (len(lon_arr), len(lat_arr))

        **Arguments:**
            *lon_arr:*
                array of longitudes (1d or 2d)
            *lat_arr:*
                array of latitudes (1d or 2d)

        **Returns:**
            *lon_arr:*
                array of longitudes (1d or 2d)
            *lat_arr:*
                array of latitudes (1d or 2d)
        """
        tmp_lon = np.array(lon_arr).squeeze()
        tmp_lat = np.array(lat_arr).squeeze()
        if np.ndim(tmp_lon) == np.ndim(tmp_lat) == 1:
            my_lon_arr = np.vstack(len(tmp_lat)*(tmp_lon,))
            my_lat_arr = np.hstack(len(tmp_lon)*(tmp_lat[:, np.newaxis],))
        else:
            my_lon_arr = tmp_lon
            my_lat_arr = tmp_lat
        return my_lon_arr, my_lat_arr

    @property
    def coordinates(self):
        """Returns coordinates

        **Returns:**
            *lon_arr:*
                longitudes
            *lat_arr:*
                latitudes
        """
        return self.lon_arr, self.lat_arr


    def get_coordinates(self):
        """Returns coordinates

        **Returns:**
            *lon_arr:*
                longitudes
            *lat_arr:*
                latitudes
        """
        return self.lon_arr, self.lat_arr

 
    def get_coordinates_geo(self):
        """Returns non-rotated geographical coordinates

        **Returns:**
            *lon_arr:*
                longitudes in geographical coordinates
            *lat_arr:*
                latitudes in geographical coordinates
        """
        if not self.rotated:
           lon_arr_geo = self.lon_arr
           lat_arr_geo = self.lat_arr
        else:
           lon_arr_geo, lat_arr_geo = self.transform().get_coordinates()  

        return lon_arr_geo, lat_arr_geo


    def get_schwedischer_Tuersteher(self):
        """Returns a Swedish bodyguard.

        **Returns:**
            *Lasse Reinström:*
                Swedish bodyguard
        """
        return "Lasse Reinström"


    def get_dimensions(self):
        """Returns the dimensions of the grid.

        **Returns:**
            *dimensions:*
                dimensions of the grid
        """
        return self.lon_arr.shape


    @property
    def pole(self):
        """Returns the coordinates of the North Pole

        **Returns:**
            *pol_lon:*
                longitude of the North Pole
            *pol_lat:*
                latitude of the North Pole
        """
        return (self.pol_lon, self.pol_lat)


    def get_boundary_as_polygon(self, do_geo=True):
        """Returns lon-lat information at the boundary.

        **Returns:**
            *lon_square:*
                list of longitudes describing the boundary-polygon
            *lat_square:*
                list of latitudes describing the boundary-polygon
        """
        xhor, yhor = self.get_coordinates()
        dimensions = xhor.shape
        xbottom = xhor[0, :]
        xright = xhor[:, dimensions[1]-1]
        xtop = xhor[dimensions[0]-1, :][::-1]
        xleft = xhor[:, 0][::-1]

        ybottom = yhor[0, :]
        yright = yhor[:, dimensions[1]-1]
        ytop = yhor[dimensions[0]-1, :][::-1]
        yleft = yhor[:, 0][::-1]

        lon_square = np.concatenate((xbottom, xright, xtop, xleft))
        lat_square = np.concatenate((ybottom, yright, ytop, yleft))

        return lon_square, lat_square


    def get_bounding_box(self):
        """Returns lon-lat information of the grid's bounding box.

        The bounding box is defined by the four corners of a rectangel enclosing
        the complete area of the grid's coordinates.

        **Returns:**
            *ll, ul, ur, lr*
               Tuples of the four corners of the bounding box
        """
        lon, lat = self.coordinates

        ll = (np.min(lon),np.min(lat))
        ul = (np.min(lon),np.max(lat))
        ur = (np.max(lon),np.max(lat))
        lr = (np.max(lon),np.min(lat))

        return (ll, ul, ur, lr)


    def get_center(self):
        """Returns lon-lat information of center grid cell.

        The coordinates of the grid's center are not neccessarily
        at (0,0).

        **Returns:**
            *lon lat*
               Centered grid cell coordinate
        """
        lon, lat = self.coordinates

        dimx = lon.shape[0]
        dimy = lon.shape[1]
 
        return (lon[dimx/2][dimy/2],lat[dimx/2][dimy/2])


    
    def get_corners(self):
        """Returns lon-lat information of the four grid corners.

        The four corners of the grid's coordinates.

        **Returns:**
            *ll, ul, ur, lr*
               Tuples of the four corners of the grid
        """
        lon, lat = self.coordinates
        
        ll = (lon[0][0],lat[0][0])
        ul = (lon[-1][0],lat[-1][0])
        ur = (lon[-1][-1],lat[-1][-1])
        lr = (lon[0][-1],lat[0][-1])

        return (ll, ul, ur, lr) 


    def get_grid_box(self, lon, lat):
        """Returns the grid box containing the given point.

        """
        lon_1d = self.lon_arr[1,:]
        lat_1d = self.lat_arr[:,1]
        box_number_x = self._ret_box_position(lon_1d, lon)
        box_number_y = self._ret_box_position(lat_1d, lat)

        return (box_number_x, box_number_y)


    def transform(self, pol_lon=None, pol_lat=None):
        """Returns a transformed Grid.

        **Attributes:**
            *pol_lon:*
                longitude of rotated North Pole
            *pol_lat:*
                latitude of rotated North Pole

        **Returns:**
            *Grid:*
                New rotated Grid.

        Written by Lars Buntemeyer
        """
        if self.rotated:
            direction = 'rot2geo'
            pol_lon = self.pol_lon
            pol_lat = self.pol_lat
        else:
            if pol_lon is None or pol_lat is None:
                pol_lon = self.pol_lon
                pol_lat = self.pol_lat
                #raise Exception('grid is not rotated, transform requires pol_lon and pol_lat')
            direction = 'geo2rot'
        lon_arr_trans, lat_arr_trans = rotated_grid_transform(
            self.lon_arr, self.lat_arr, pol_lon, pol_lat,
            direction=direction)
        if self.rotated:
            return Grid(lon_arr_trans, lat_arr_trans)
        else:
            return Grid(lon_arr_trans, lat_arr_trans, pol_lon, pol_lat)


    @staticmethod
    def _ret_box_position(coord_arr, coord):
        """
        """
        min_distance = 360.0
        box_number = 0
        for position, coordinate in enumerate(coord_arr, start=1):
            distance = abs(coordinate - coord)
            if distance < min_distance:
                min_distance = distance
                box_number = position
            elif distance > min_distance:
                break

        return box_number



class RotGrid(Grid):
    """This class contains rotated or non-rotated grid-information.

    We are looking at the same points in different coordinate-systems:

    * Variables ending with _geo are describing the points in the
      geographical non-rotated coordinates system you find on a globe.

    * Variables ending with _rot are describing the same points in the
      rotated coordinates system.

    **Attributes:**
        *lon_arr_rot:*
            longitudes in rotated coordinate system (2d-array)
        *lat_arr_rot:*
            latitudes in rotated coordinate system (2d-array)
        *pol_lon_geo:*
            longitude of rotated North Pole
        *pol_lat_geo:*
            latitude of rotated North Pole
        *lon_arr_geo:*
            longitudes in geographical system (2d-array)
        *lat_arr_geo:*
            latitudes in geographical system (2d-array)

    .. note::
        RotGrid class is deprecated. Please use Grid in the future.
    """


    # Methods
    def __init__(self, lon_arr, lat_arr, pol_lon, pol_lat):
        """Setting lon/lat-array and rotated North Pole
        
        **Arguments:**
            *lon_arr:*
                longitudes in rotated coordinate system (2d-array)
            *lat_arr:*
                latitudes in rotated coordinate system (2d-array)
            *pol_lon:*
                longitude of rotated North Pole
            *pol_lat:*
                latitude of rotated North Pole
        """
        logging.warning('RotGrid class will not be supported in the future. Please create '+
              'a Grid instance instead with a rotated pole if required.')
        Grid.__init__(self, lon_arr, lat_arr, pol_lon, pol_lat)


    @property
    def lon_arr_rot(self):
        return self.lon_arr

    @property
    def lat_arr_rot(self):
        return self.lat_arr

    @property
    def pol_lon_geo(self):
        return self.pol_lon

    @property
    def pol_lat_geo(self):
        return self.pol_lat

    def get_coordinates_rot(self):
        """Returns rotated coordinates

        Subroutine is here for downwards compatibility.
        
        **Returns:**
            *lon_rot_arr:*
                longitudes in rotated coordinates
            *lat_rot_arr:*
                latitudes in rotated coordinates
        """
        return self.get_coordinates()


    def get_rot_pole(self):
        """Returns the coordinates of the rotated North Pole

        Subroutine is here for downwards compatibility.


        **Returns:**
            *pol_lon:*
                longitude of the rotated North Pole
            *pol_lat:*
                latitude of the rotated North Pole
        """
        return self.pole




def rotated_coord_transform(lon, lat, np_lon, np_lat,
                            direction='rot2geo'):
    """Transforms a coordinate into a rotated grid coordinate and vice versa.

    The coordinates have to given in degree and will be returned in degree.

    **Arguments:**
        *lon:*
            Longitude coordinate.
        *lat:*
            Latitude coordinate.
        *np_lon:*
            Longitude coordinate of the rotated pole.
        *np_lat:*
            Latitude coordinate of the rotated pole.
        *direction:*
            Direction of the rotation.
            Options are: 'rot2geo' (default) for a transformation to regular
            coordinates from rotated. 'geo2rot' transforms regular coordinates
            to rotated.

    **Returns:**
        *lon_new:*
            New longitude coordinate.
        *lat_new:*
            New latitude coordinate.

    Written by Kevin Sieck
    """

    # Convert degrees to radians
    lon = (lon * math.pi) / 180.
    lat = (lat * math.pi) / 180.

#    SP_lon = SP_coor(1)
#    SP_lat = SP_coor(2)

    theta = 90. - np_lat # Rotation around y-axis
    phi = np_lon + 180.  # Rotation around z-axis

    # Convert degrees to radians
    phi = (phi * math.pi) / 180.
    theta = (theta * math.pi) / 180.

    # Convert from spherical to cartesian coordinates
    x = math.cos(lon) * math.cos(lat)
    y = math.sin(lon) * math.cos(lat)
    z = math.sin(lat)

    # Regular -> Rotated
    if direction == 'geo2rot':

        x_new = (math.cos(theta) * math.cos(phi) * x +
                 math.cos(theta) * math.sin(phi) * y +
                 math.sin(theta) * z)
        y_new = (- math.sin(phi) * x +
                   math.cos(phi) * y)
        z_new = (- math.sin(theta) * math.cos(phi) * x -
                   math.sin(theta) * math.sin(phi) * y +
                   math.cos(theta) * z)

    # Rotated -> Regular
    elif direction == 'rot2geo':
    
        phi = - phi
        theta = - theta
    
        x_new = (math.cos(theta) * math.cos(phi) * x +
                 math.sin(phi) * y +
                 math.sin(theta) * math.cos(phi) * z)
        y_new = (- math.cos(theta) * math.sin(phi) * x +
                   math.cos(phi) * y -
                   math.sin(theta) * math.sin(phi) * z)
        z_new = (- math.sin(theta) * x +
                   math.cos(theta) * z)

    # Convert cartesian back to spherical coordinates
    lon_new = math.atan2(y_new, x_new)
    lat_new = math.asin(z_new)

    # Convert radians back to degrees
    lon_new = (lon_new * 180.) / math.pi
    lat_new = (lat_new * 180.) / math.pi;

    return (lon_new, lat_new)


def rotated_grid_transform(lon_arr, lat_arr, np_lon, np_lat,
                           direction='rot2geo'):
    """Transforms a grid into a rotated grid and vice versa.

    The grid coordinates have to given in degree and will be returned in degree.

    **Arguments:**
        *lon_arr:*
            Array with longitude coordinates (at least 2D).
        *lat_arr:*
            Array with latitude coordinates (at least 2D).
        *np_lon:*
            Longitude coordinate of the rotated pole.
        *np_lat:*
            Latitude coordinate of the rotated pole.
        *direction:*
            Direction of the rotation.
            Options are: 'rot2geo' (default) for a transformation to regular
            coordinates from rotated. 'geo2rot' transforms regular coordinates
            to rotated.

    **Returns:**
        *lon_arr_new:*
            New longitude coordinates.
        *lat_arr_new:*
            New latitude coordinates.

    Written by Kevin Sieck
    """
    dimensions = lon_arr.shape
    lon_arr_new = np.zeros(dimensions)
    lat_arr_new = np.zeros(dimensions)
    for i in range(dimensions[1]):       # x-direction
        for j in range(dimensions[0]):   # y-direction
            lon = lon_arr[j, i]
            lat = lat_arr[j, i]
            lon_new, lat_new = rotated_coord_transform(
                lon, lat, np_lon, np_lat, direction=direction)
            lon_arr_new[j, i] = lon_new
            lat_arr_new[j, i] = lat_new

    return (lon_arr_new, lat_arr_new)



def get_real_coord(phis, rlas, polphi, pollam):
    '''Returns the regular lat/lon coordinate of a rotated coordinate.

    This definition was taken from the REMO model and translated
    into python code.

    **Arguments:**
        *phis:*
            Latitude coordinate of the rotated grid.
        *rlas:*
            Longitude coordinate of the rotated grid.
        *polphi:*
            Latitude coordinate of the rotated pole.
        *pollam:*
            Longitude coordinate of the rotated pole.

    **Returns:**
        *(phstoph, rlstorl):*
            Tuple of the regular coordinates (lat, lon)

    Written by Kevin Sieck

    Last changes 31.10.2010
    '''
    rpi18 = 57.2957795
    pir18 = 0.0174532925

    sinpolp = math.sin(pir18*polphi)
    cospolp = math.cos(pir18*polphi)

    sinpoll = math.sin(pir18*pollam)
    cospoll = math.cos(pir18*pollam)

    sinphis = math.sin(pir18*phis)
    cosphis = math.cos(pir18*phis)

    if rlas > 180.0:
        rlas = rlas - 360.0

    sinrlas = math.sin(pir18*rlas)
    cosrlas = math.cos(pir18*rlas)

    # compute latitude coordinate
    arg = cospolp*cosphis*cosrlas + sinpolp*sinphis
    phstoph = rpi18*math.asin(arg)

    # compute longitude coordinate
    arg1 = sinpoll*(- sinpolp*cosrlas*cosphis +
                      cospolp*sinphis) - cospoll*sinrlas*cosphis
    arg2 = cospoll*(- sinpolp*cosrlas*cosphis +
                      cospolp*sinphis) + sinpoll*sinrlas*cosphis

    if arg2 == 0.0:
        arg2 = pow(10, -20)

    rlstorl = rpi18*math.atan2(arg1, arg2)

    return(phstoph, rlstorl)




def from_cdo_griddes(griddes):
    """Returns a Grid instance from reading a cdo griddes file.

    **Arguments:**
        *griddes:*
            Textfile with a grid description (cdo style).

    **Returns:**
        *grid:*
            Grid instance.

    Written by Kevin Sieck
    """

    with open(griddes) as grid_file:
        grid_file_lines = grid_file.readlines()

    grid_dic = {}

    for line in grid_file_lines:
        words = line.split()
        if words[0] == '#':
            continue
        else:
            length = len(words)
            if length == 3:
                grid_dic[words[0]] = words[2]
            else:
                value_string = ' '.join(words[2:length-1])
                grid_dic[words[0]] = value_string

    if grid_dic['gridtype'] != 'lonlat':
        print(('Gridtype {0} not supported'.format(grid_dic['gridtype'])))
        return ''

    lon = np.zeros(int(grid_dic['xsize']))
    lat = np.zeros(int(grid_dic['ysize']))

    for i in range(len(lon)):
        lon[i] = float(grid_dic['xfirst']) + i * float(grid_dic['xinc'])
    for j in range(len(lat)):
        lat[j] = float(grid_dic['yfirst']) + j * float(grid_dic['yinc'])

    if grid_dic['xname'] == 'rlon':
        pol_lon = float(grid_dic['xnpole'])
        pol_lat = float(grid_dic['ynpole'])
        grid = RotGrid(lon, lat, pol_lon, pol_lat)
    else:
        grid = Grid(lon, lat)

    return grid
