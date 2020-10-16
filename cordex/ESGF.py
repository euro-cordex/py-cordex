# -*- coding: utf-8 -*-
# flake8: noqa
"""ESGF conventions module

This module defines common ESGF path and filename conventions.

Example:

    To get a list of available implementations, you can call, e.g.,::

        from cordex import ESGF
        print(ESGF.conventions())

Example:

    To get an instance of the :class:`ESGFFileSelection` you can use
    the main interface function::

        example of interface function here.

"""

import os
import copy
import logging
import datetime as dt
import pandas as pd

from . import conventions as conv


cordex_path_list = ['product','CORDEX_domain','institute_id','driving_model_id', \
                    'experiment_id', 'ensemble_member', 'model_id', 'rcm_version_id'  , \
                    'frequency', 'variable', 'version']
cordex_conv_str  = '{variable}_{CORDEX_domain}_{driving_model_id}_{experiment_id}_' \
                   '{ensemble_member}_{model_id}_{rcm_version_id}_{frequency}_' \
                   '{startdate}-{enddate}.{suffix}'


UNIQUE = ['product', 'CORDEX_domain', 'institute_id', 'driving_model_id', 'experimentd_id',
          'ensemble_member', 'model_id', 'rcm_version_id', 'frequency', 'variable', 'version', 'modeling_realm',
          'mip_table', 'variable']

date_fmts = {12:'%Y%m%d%H%M', 10:'%Y%m%d%H', 8:'%Y%m%d', 6:'%Y%m', 4:'%Y'}



def parse_date(date_str):
    return dt.datetime.strptime(date_str, date_fmts[len(date_str)])

def format_date(date, freq):
    return date.strftime(date_fmts[freq])



class ESGFFileNameConvention(conv.FileNameConvention):

    def __init__(self, *args, **kwargs):
        conv.FileNameConvention.__init__(self, *args, **kwargs)

    #def parse_attrs(self, attrs):
    #    print(attrs)
    #    attrs['startdate'] = parse_date(attrs['startdate'], attrs['frequency'])
    #    attrs['enddate']   = parse_date(attrs['enddate'], attrs['frequency'])
    #    return attrs

    #def format_attrs(self, attrs, any_str):
    #    if isinstance(attrs['startdate'], dt.datetime):
    #        attrs['startdate'] = format_date(attrs['startdate'], attrs['frequency'])
    #    if isinstance(attrs['enddate'], dt.datetime):
    #        attrs['enddate']   = format_date(attrs['enddate'], attrs['frequency'])
    #    return attrs


class CORDEX(conv.FileConvention):
    """Implements CORDEX path and filename conventions.
    """
    name      = 'CORDEX'

    def __init__(self, root=None):
        path_conv      = conv.FilePathConvention(cordex_path_list, root=root)
        filename_conv  = ESGFFileNameConvention(cordex_conv_str)
        conv.FileConvention.__init__(self, path_conv, filename_conv)



class CMIP5():
    """Implements CMIP5 path and filename conventions.
    """
    name      = 'CMIP5'
    cmip5_path_list  = ['product','institute','model', \
                        'experiment', 'frequency', 'modeling_realm', 'mip_table', 'ensemble_member'  , \
                        'version', 'variable']
    cmip5_conv_dyn   = '{variable}_{mip_table}_{model}_{experiment}_' \
                       '{ensemble_member}_{startdate}-{enddate}.{suffix}'
    cmip5_conv_fx    = '{variable}_{mip_table}_{model}_{experiment}_' \
                       '{ensemble_member}.{suffix}'
    fx_vars = ['orog', 'sftlf', 'sftof', 'areacella', 'volcello', 'deptho', 'areacello']

    def __init__(self, root=None):
        path_conv      = conv.FilePathConvention(self.cmip5_path_list, root=root)
        filename_dyn   = ESGFFileNameConvention(self.cmip5_conv_dyn)
        filename_fx    = ESGFFileNameConvention(self.cmip5_conv_fx)
        self.conv_dyn = conv.FileConvention(path_conv, filename_dyn)
        self.conv_fx  = conv.FileConvention(path_conv, filename_fx)

    def parse(self, file):
        """Parses a file including path and filename and returns attributes.
        """
        path_attrs     = self.conv_dyn.path_conv.parse(os.path.dirname(file))
        if path_attrs['variable'] in self.fx_vars:
            attrs = self.conv_fx.parse(file)
        else:
            attrs = self.conv_dyn.parse(file)
        return attrs

    def pattern(self, **kwargs):
        if 'variable' in kwargs:
            if kwargs['variable'] in self.fx_vars:
                return self.conv_fx.pattern(**kwargs)
            else:
                return self.conv_dyn.pattern(**kwargs)



class _ConventionFactory(object):
    """Factory class for creating a NamingConvention instance.
    """

#    @staticmethod
#    def create_convention(name):
#        path_conv     = conv.FilePathConvention(ESGF_CONVS[name]['path'])
#        filename_con  = conv.FileNameConvention(ESGF_CONVS[name]['file'])
#        return conv.FileConvention(path_conv, filename_conv)

    @staticmethod
    def conventions():
#        return [cls.create_convention(name) for name in ESGF_CONVS]
        return [CORDEX, CMIP5]

    @classmethod
    def names(cls):
        list_of_names = []
        for conv in cls.conventions():
           list_of_names.append(conv.name)
        return list_of_names

    @classmethod
    def get_convention(cls, name):
        convention = None
        for conv in cls.conventions():
           if name == conv.name:
             convention = conv
        if convention is None:
           logging.error('Unknown convention name: '+name)
           logging.info('Known convention names: '+str(cls.names()))
           raise Exception('Unknown convention name: '+name)
        else:
           return convention



class ESGFFileSelection(conv.FileSelection):
    """Holds a Pandas Dataframe object.

    The :class:`ESGFFileSelection` holds and manages a Pandas
    Dataframe instance. It defines some common methods to work
    with ESGF netcdf files.
    """
    def __init__(self, *args, **kwargs):
        conv.FileSelection.__init__(self, *args, **kwargs)

    def __str__(self):
        text = ''
        text += super().__str__() + '\n'
        for col in UNIQUE:
            if col in self.df:
                text += "{:<30}  :   {}\n".format(col, self.df[col].unique())
        text += "{:<30}  :   {} to {}\n".format('time range', self.timerange[0], self.timerange[1])
        if self.unique:
            text += '\nESGF File Selection is unique.\n'
        else:
            text += '\nESGF File Selection is not unique.\n'
        return text

    def subset(self, **kwargs):
        """Creates a subset by filtering attributes.
        """
        return ESGFFileSelection(super().subset(**kwargs).df)

    def to_datetime(self):
        """Converts the date columns to datetime objects.

        The date columns (startdate, enddate) are converted to datetime
        objects depending on the lenght of the date string.

        Returns:
             :class:`ESGFFileSelection`: selection converted date columns.
        """
        df = self.df
        for index, row in df.iterrows():
            row['startdate'] = parse_date(row['startdate'])
            row['enddate']   = parse_date(row['enddate'])
        df.sort_values(by='startdate', inplace=True)
        return ESGFFileSelection(df)

    def select_timerange(self, time_range):
        """Returns a selected timerange.

        Args:
            time_range (tuple): Tuple that contains a startdate
                and enddate in datetime format.

        Returns:
             :class:`ESGFFileSelection`: selection within time range.

        """
        logging.debug('selecting time range: {} tp {}'.format(time_range[0], time_range[1]))
        df = self.df[((self.df['startdate'] >= time_range[0]) & (self.df['enddate'] <= time_range[1])) |
                ((self.df['startdate'] <= time_range[0]) & (self.df['enddate'] >= time_range[0])) |
                ((self.df['startdate'] <= time_range[1]) & (self.df['enddate'] >= time_range[1]))]
        return ESGFFileSelection(df)

    @property
    def timerange(self):
        if 'startdate' in self.df and 'enddate' in self.df:
            return (min(self.df['startdate']), max(self.df['enddate']))
        else:
            return (None, None)

    @property
    def unique(self):
        """Determines if the file selection is unique.

        The file selection is unique, if all columns except the start
        and end dates have single unique values. If that is the case,
        the file selection holds a time series of single unique attributes.

        This is useful to check if you have selected a time series of a
        single variable that might be processed or plotted.

        Returns:
            True if selection is unique, false otherwise.

        """
        for column, data in self.df.items():
            column_in_unique = column in UNIQUE
            data_unique      = len(data.unique())==1
            if column_in_unique and not data_unique:
                return False
        return True


def select_files(project_id, filter={}, root=None, **kwargs):
    """Creates a file selection containing attributes.
    """
    convention = get_convention(project_id, root=root)
    return conv.select_files(convention, filter, root, **kwargs)


def file_selection_from_scratch():
    pass

def file_selection_from_csv(filename):
    return ESGFFileSelection(pd.from_csv(filename))

def get_selection(convention_id, filter={}, root=None, **kwargs):
    """Top level function to create a :class:`ESGFFileSelection` instance.

    This function creates a :class:`FileSelection` instance
    using a file naming convention of type :class:``FileConvention`.

    Args:
        convention_id (str): The name of the convention.
        filter (dict): Defines attributes to filer the search.
        root (str): The root directory where the convention holds.

    Returns:
        :class:`ESGFFileSelection` object.

    """
    convention = get_convention(convention_id, root=root)
    files      = conv.select_files(convention, filter, root, **kwargs)
    df         = conv.make_df(convention, files)
    return ESGFFileSelection(df)


def conventions():
    """Lists available ESGF conventions.

    Returns:
        List of available ESGF convention ids.
    """
    return _ConventionFactory.names()


def get_convention(name, root=None):
    """Returns a ESGS convention instance.

    Args:
        name (str): The convention id.

    Returns:
        :class:`ESGFFileConvention` object.
    """
    return _ConventionFactory.get_convention(name)(root=root)




def iterdict(d):
    for k,v in d.items():
        if isinstance(v, dict):
            iterdict(v)
        else:
            print (k,":",v)



        #for item in reversed(path_conv[:-1]):
        #    dict = {item:dict} 
        #for path in self.pathes:
        #    path_list = path.split(os.sep)
        #    for value in reversed(path_list[-nitem:]):


