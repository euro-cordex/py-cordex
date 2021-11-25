import pandas as pd

from ._resources import (
    read_cordex_domain_tables,
    ecmwf_tables,
    fetch_cordex_cmor_table,
    cmor_table_version,
)

__cmor_table_version__ = cmor_table_version


class read_cls:
    def __init__(self, reader):
        self.reader = reader

    @property
    def tables(self):
        return self.reader()

    @property
    def table(self):
        return pd.concat(self.tables.values())

    # def __getattr__(self, table):
    #    return self.tables[table]


domains = read_cls(read_cordex_domain_tables)

ecmwf = read_cls(ecmwf_tables)


def cordex_cmor_table(table, table_dir=None):
    """fetch a cordex cmor table

    If required, the table will be download from github.
    The tables are experimental right now and only used
    for development purposes.

    Parameters
    ----------
    table: str
        Name of the cordex table without the CORDEX_ prefix
        and the .json suffix.
    table_dir: str
        Path to table directory. If `None`, the tables will
        be downloaded from the github repo if neccessary.

    Returns
    -------
    filename : str
        Filepath to the Cordex cmor table.
    """
    return fetch_cordex_cmor_table(table)
