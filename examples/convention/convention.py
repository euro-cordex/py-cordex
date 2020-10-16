
from cordex import conventions as conv

# create your own filename convention string and list
filename_conv_str  = 'my_convention_{variable}_{model_id}_{domain_id}.nc'
path_conv_list     = ['model_id','variable']

# create conventions for filename and path
filename_conv = conv.FileNameConvention(filename_conv_str)
path_conv     = conv.FilePathConvention(path_conv_list)


# now define your attributes to fill the templates.
root = '/my_root'
attributes = {'model_id'        : 'GERICS-REMO2015',
              'variable'        : 'pr',
              'domain_id'       : 'EUR-11'}

# create filename and path
filename = filename_conv.pattern(**attributes)
path     = path_conv.pattern(root, **attributes)

# create combined file convention
file_conv = conv.FileConvention(path_conv, filename_conv)

# create full filename with path
file = file_conv.pattern(root, **attributes)
