# -*- coding: utf-8 -*-
# flake8: noqa

from . import dbase as db


CF_URL = 'http://cfconventions.org/Data/cf-standard-names/{}/src/cf-standard-name-table.xml'


def get_cf_table(version='70', url=CF_URL):
    return db.getxml(CF_URL.format(version))



#char rotated_pole ;
#        rotated_pole:grid_mapping_name = "rotated_latitude_longitude" ;
#        rotated_pole:grid_north_pole_latitude = 39.25f ;
#        rotated_pole:grid_north_pole_longitude = -162.f ;
#double rlon(rlon) ;
#        rlon:axis = "X" ;
#        rlon:standard_name = "grid_longitude" ;
#        rlon:long_name = "longitude in rotated pole grid" ;
#        rlon:units = "degrees" ;
#double lon(rlat, rlon) ;
#        lon:standard_name = "longitude" ;
#        lon:long_name = "longitude" ;
#        lon:units = "degrees_east" ;
#double rlat(rlat) ;
#        rlat:axis = "Y" ;
#        rlat:standard_name = "grid_latitude" ;
#        rlat:long_name = "latitude in rotated pole grid" ;
#        rlat:units = "degrees" ;
#double lat(rlat, rlon) ;
#        lat:standard_name = "latitude" ;
#        lat:long_name = "latitude" ;
#        lat:units = "degrees_north" ;


coords = \
        {
   'rlon'         : {'axis': 'X',
                     'standard_name': 'grid_longitude',
                     'long_name'    : 'longitude in rotated pole grid',
                     'units'        : 'degrees'},
   'rlat'         : {'axis': 'Y',
                     'standard_name': 'grid_latitude',
                     'long_name'    : 'latitude in rotated pole grid',
                     'units'        : 'degrees'},
   'lon'          : {
                     'standard_name': 'longitude',
                     'long_name'    : 'longitude',
                     'units'        : 'degrees_east'},
   'lat'          : {
                     'standard_name': 'latitude',
                     'long_name'    : 'latitude',
                     'units'        : 'degrees_north'},
   }

mapping = \
     { 'grid_mapping_name' : 'rotated_latitude_longitude',
                     'grid_north_pole_latitude' : 90.,
                     'grid_north_pole_longitude': 180.,
                     'north_pole_grid_longitude': 0.}


DEFAULT_MAPPING_NCVAR = 'rotated_latitude_longitude'
