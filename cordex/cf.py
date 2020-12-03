# -*- coding: utf-8 -*-
# flake8: noqa



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
