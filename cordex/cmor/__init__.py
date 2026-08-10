from .cmor import Cmorizer, cmorize_variable, prepare_variable
from .config import options, set_options
from .utils import (
    mid_of_month,
    mid_of_season,
    month_bounds,
    season,
    season_bounds,
    to_cftime,
)

# def fetch_basic_tables():
#     """fetch basic cmor tables"""
#     tables.cmip6_cmor_table("CORDEX_coordinate.json")
#     tables.cmip6_cmor_table("CORDEX_grids.json")
#     tables.cmip6_cmor_table("CORDEX_formula_terms.json")
#     tables.cordex_cmor_table("CORDEX_CV.json")


# fetch_basic_tables()


__all__ = [
    "Cmorizer",
    "cmorize_variable",
    "mid_of_month",
    "mid_of_season",
    "month_bounds",
    "options",
    "prepare_variable",
    "season",
    "season_bounds",
    "set_options",
    "to_cftime",
]
