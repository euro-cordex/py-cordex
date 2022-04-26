import pytest

import cordex as cx


@pytest.mark.parametrize(
    "table",
    [
        "CMIP6_Amon",
        "CMIP6_day",
        "CMIP6_E1hr",
        "CMIP6_E3hr",
        "CMIP6_coordinate",
        "CMIP6_fx",
        "CMIP6_grids",
    ],
)
def test_download_cmip6_cmor_tables(table):
    cx.tables.cmip6_cmor_table(table)


@pytest.mark.parametrize("table", ["CORDEX_CV", "CORDEX_remo_example"])
def test_download_cordex_cmor_tables(table):
    cx.tables.cordex_cmor_table(table)
