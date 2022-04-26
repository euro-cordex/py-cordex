from .. import tables
from .utils import (
    mid_of_month,
    mid_of_season,
    month_bounds,
    season,
    season_bounds,
    to_cftime,
)


def fetch_basic_tables():
    """fetch basic cmor tables"""
    tables.cmip6_cmor_table("CMIP6_coordinate.json")
    tables.cmip6_cmor_table("CMIP6_grids.json")
    tables.cmip6_cmor_table("CMIP6_formula_terms.json")
    tables.cordex_cmor_table("CORDEX_CV.json")


fetch_basic_tables()


__all__ = [
    "mid_of_month",
    "mid_of_season",
    "month_bounds",
    "season",
    "season_bounds",
    "to_cftime",
]
