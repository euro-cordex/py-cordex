from .utils import (
    to_cftime,
    season,
    season_bounds,
    mid_of_season,
    mid_of_month,
    month_bounds,
    to_cftime,
)


from .. import tables


def fetch_basic_tables():
    """fetch basic cmor tables"""
    tables.cmip6_cmor_table('CMIP6_coordinate.json')
    tables.cmip6_cmor_table('CMIP6_grids.json')
    tables.cmip6_cmor_table('CMIP6_formula_terms.json')
    tables.cordex_cmor_table('CORDEX_CV.json')
    
fetch_basic_tables()