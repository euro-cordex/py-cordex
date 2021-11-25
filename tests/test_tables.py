import cordex as cx
import pytest


@pytest.mark.parametrize(
    "table", ["Amon", "day", "1hr", "3hr", "CV", "coordinate", "fx", "grids"]
)
def test_download_tables(table):
    cx.cordex_cmor_table(table)
