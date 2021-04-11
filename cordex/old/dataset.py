

from netCDF4 import Dataset, MFDataset, date2index, num2date, date2num
import xarray as xr


class NC4Dataset():

    def __init__(self, file_list=None, ds=None, time_axis='time', **kwargs):
        self.time_axis = time_axis
        if file_list:
            self.ds = Dataset(file_list, **kwargs)
        else:
            self.ds = ds

    def __str__(self):
        return str(self.ds)

    @property
    def calendar(self):
        return self.ds.variables[self.time_axis].calendar

    @property
    def units(self):
        return self.ds.variables[self.time_axis].units

    @property
    def variables(self):
        return self.ds.variables

    def __getattr__(self, item):
        return getattr(self.ds, item)

    def ncattrs(self):
        return self.ds.ncattrs()

    def ncattrs_dict(self, varname=None):
        if varname:
            return {attr:self.ds.variables[varname].getncattr(attr)
                    for attr in self.ds.variables[varname].ncattrs()}
        else:
            return {attr:self.ds.getncattr(attr)
                    for attr in self.ds.ncattrs()}

    def get_index_by_date(self, dates):
        return date2index(dates,
                self.ds.variables[self.time_axis], select='exact')

    def get_date_by_num(self, nums):
        return num2date(nums, self.units, calendar=self.calendar)

    def get_num_by_date(self, dates):
        return date2num(dates, self.units, calendar=self.calendar)

    def get_num_by_index(self, index):
        return self.variables[self.time_axis][index]

    def get_date_by_index(self, index):
        num = self.ds.variables[self.time_axis][index]
        return self.get_date_by_num(num)

    def get_timestep(self, date, varname):
        ix = self.get_index_by_date(date)
        return self.ds.variables[varname][ix]

    def get_variable(self, varname):
        return self.ds.variables[varname]

    def dynamic(self, varname):
        return self.ds.variables[varname].dimensions[0] == self.time_axis

    def data_by_date(self, variable, date):
        return self.timestep(variable, self.get_index_by_date(date))

    def static(self, varname):
        return not self.dynamic(varname)

    def timestep(self, varname, timestep=0):
        if self.dynamic(varname):
            return self.variables[varname][timestep]
        else:
            return self.variables[varname][:]

    def threeD(self, varname):
        return len(self.timestep(varname).shape) == 3

    def twoD(self, varname):
        return not self.threeD(varname)

    def nlev(self, varname=None):
        if self.threeD(varname):
            return self.ds.variables[varname].shape[1:-2][0]
        else:
            return 1

    def to_netcdf(self, varname=None, timestep=None, destination=None):
        ds = copy_dataset(self.ds, varname=varname, timestep=timestep, destination=destination)
        filepath = ds.filepath()
        ds.close()
        return filepath




class NC4MFDataset(NC4Dataset):

    def __init__(self, file_list, time_axis='time', **kwargs):
        self.time_axis = time_axis
        self.ds = MFDataset(file_list, **kwargs)

    def ncattrs_dict(self, varname=None):
        if varname:
            return {attr:getattr(self.ds.variables[varname], attr)
                    for attr in self.ds.variables[varname].ncattrs()}
        else:
            return {attr:getattr(self.ds, attr)
                    for attr in self.ds.ncattrs()}

    def getncattr(self, attr):
        return getattr(self.ds, attr)



class XRDataset():

    def __init__(self, file_list, time_axis='time'):
        self.time_axis = time_axis
        self.ds = xr.open_dataset(file_list, decode_times=False, decode_cf=False,
                decode_coords=False, chunks={'time':1, 'lat':1000, 'lon':1000})


class XRMFDataset():

    def __init__(self, file_list, time_axis='time'):
        self.time_axis = time_axis
        self.ds = xr.open_mfdataset(file_list)


def copy_dataset(src, varname=None, timestep=None, destination=None):
    if varname is None:
        variables = src.variables
    else:
        variables = {varname:src.variables[varname]}
    if destination is None:
        destination = 'copy.nc'
    dst = Dataset(destination,'w')
    # copy attributes
    for name in src.ncattrs():
        dst.setncattr(name, getattr(src, name))
    # copy dimensions
    for name, dimension in src.dimensions.items():
        if timestep and name=='time':
            length = 1
        else:
            length = (len(dimension) if not dimension.isunlimited() else None)
        #dst.createDimension( name, (len(dimension) if not dimension.isunlimited() else None))
        dst.createDimension( name, length)
        # copy all file data except for the excluded
    for name, variable in variables.items():
        if hasattr(variable, '_FillValue'):
            fill_value = getattr(variable, '_FillValue')
        else:
            fill_value = None
        var = dst.createVariable(name, variable.dtype, variable.dimensions, fill_value=fill_value)
        for attr in variable.ncattrs():
            if attr != '_FillValue':
                dst.variables[name].setncattr(attr, getattr(variable, attr))
        if variable.shape:
            if timestep is None or 'time' not in variable.dimensions:
                dst.variables[name][:] = src.variables[name][:]
            else:
                dst.variables[name][0] = src.variables[name][timestep]
    return dst
