from cordex import ESGF as esgf

root       = '/my_root'

# define attributes
attributes   = {'institute_id'    : 'GERICS',
                'product'         : 'output',
                'model_id'        : 'GERICS-REMO2015',
                'experiment_id'   : 'evaluation',
                'driving_model_id': 'ECMWF-ERAINT',
                'variable'        : 'pr',
                'rcm_version_id'  : 'v1',
                'date'            : 'v20200221',
                'frequency'       : 'day',
                'CORDEX_domain'   : 'EUR-11',
                'suffix'          : 'nc',
                'ensemble'        : 'r1i1p1'}

# we use the CORDEX convention as example
convention = esgf.CORDEX()
# print the convention patterns
print(convention.path_conv.conv_str)
print(convention.filename_conv.conv_str)
# only filename
filename = convention.filename(**attributes, startdate='20010101', enddate='20010131')
# only path
path     = convention.path(**attributes, startdate='20010101', enddate='20010131')
# only filename with path
file     = convention.pattern(root, **attributes, startdate='20010101', enddate='20010131')
