"""This module defines the csv tables for cordex.
"""

import pandas as pd
import pkg_resources

from . import domains as dm
from . import data_request as dr

#domain_tables_external = {'cordex-core': 'https://raw.githubusercontent.com/euro-cordex/tables/master/domains/cordex-core.csv',
#        'cordex': 'https://raw.githubusercontent.com/euro-cordex/tables/master/domains/cordex.csv',
#        'cordex-high-res': 'https://raw.githubusercontent.com/euro-cordex/tables/master/domains/cordex-fps.csv',
#        'cordex-fps': 'https://raw.githubusercontent.com/euro-cordex/tables/master/domains/cordex-high-res.csv'}

#data_request_tables = {'cmip5': 'https://raw.githubusercontent.com/euro-cordex/tables/master/data-request/cordex-cmip5.csv'}


def read_table(table, **kwargs):
    """reads a csv table from an external resource.
    """
    csv_file = table
    return pd.read_csv(csv_file, **kwargs)

def read_tables(csv_dict, **kwargs):
    """reads all tables from an external resource.
    """
    tables = {}
    for set_name, table in csv_dict.items():
        tables[set_name] = read_table(table, **kwargs)
    return tables

def read_resource_table(resource, csv, **kwargs):
    """reads a csv table from the package resource.
    """
    csv_file = pkg_resources.resource_stream(resource, csv)
    return pd.read_csv(csv_file, **kwargs)

def read_resource_tables(resource, csv_dict, **kwargs):
    """reads all csv tables from the package resource.
    """
    tables = {}
    for set_name, table in csv_dict.items():
        tables[set_name] = read_resource_table(resource, table, **kwargs)
    return tables



domains   = read_resource_tables('cordex.tables.domains', dm.tables, index_col='short_name')
variables = read_resource_tables('cordex.tables.data_request', dr.tables, index_col='variable_id', converters = {'frequency': (lambda x: x.split(",")) })
