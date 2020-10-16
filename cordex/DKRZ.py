

from . import conventions as conv


# era interim: EIsf00_freq_date_code
# era5: /pool/data/ERA5/ml00_1H/1979/E5ml00_1H_1979-01-01_129


class ERA5(conv.FileConvention):
    """Implements ERA5 path and filename convention at DKRZ.
    """
    name      = 'ERA5'
    era5_path_list  = ['subdir','year']
    era5_conv_str   = '{id:.2}{leveltype:.2}{hour:02d}_{frequency}_{date}_{code:.3}'


    def __init__(self, root=None):
        path_conv      = conv.FilePathConvention(self.era5_path_list, root=root)
        filename_conv  = conv.FileNameConvention(self.era5_conv_str)
        conv.FileConvention.__init__(self, path_conv, filename_conv)


class ECMWF(conv.FileConvention):
    """Implements ECMWF path and filename convention at DKRZ.
    """
    name      = 'ECMWF'
    era5_path_list  = ['subdir','year']
    era5_conv_str   = '{id:.2}{leveltype:.2}{hour:02d}_{frequency}_{date}_{code:.3}'


    def __init__(self, root=None):
        path_conv      = conv.FilePathConvention(self.era5_path_list, root=root)
        filename_conv  = conv.FileNameConvention(self.era5_conv_str)
        conv.FileConvention.__init__(self, path_conv, filename_conv)
