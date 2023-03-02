from warnings import warn

import pandas as pd

from ._resources import (
    cmor_tables_inpath,
    ecmwf_tables,
    fetch_cmip6_cmor_table,
    fetch_cordex_cmor_table,
    read_cordex_domain_tables,
)

# __cmor_table_version__ = cmor_table_version
__all__ = [
    "cmor_tables_inpath",
    "ecmwf_tables",
    "fetch_cmpi6_cmor_table",
    "fetch_cordex_cmor_table",
    "read_cordex_domain_tables",
    "cmor_tables_inpath",
]


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
    """fetch a experimental cordex cmor table

    If required, the table will be download from github.
    The tables are experimental right now and only used
    for development purposes.

    Parameters
    ----------
    table: str
        Name of the cordex table.

    Returns
    -------
    filename : str
        Filepath to the cordex cmor table.
    """
    # make sure these are fetched because cmor
    # requires them silently
    fetch_cordex_cmor_table("CORDEX_coordinate")
    fetch_cordex_cmor_table("CORDEX_grids")
    fetch_cordex_cmor_table("CORDEX_formula_terms")
    fetch_cordex_cmor_table("CORDEX_CV")

    return fetch_cordex_cmor_table(table)


def cmip6_cmor_table(table):
    """fetch a cmip6 cmor table

    If required, the table will be download from github.
    The tables are experimental right now and only used
    for development purposes.

    Parameters
    ----------
    table: str
        Name of the cmip6 table.

    Returns
    -------
    filename : str
        Filepath to the cmip6 cmor table.
    """
    warn(
        "CMIP6 cmor table fetching is deprecated and will be removed in the future. Please use cordex_cmor_table instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fetch_cmip6_cmor_table(table)
