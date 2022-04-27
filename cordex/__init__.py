# -*- coding: utf-8 -*-
# from pkg_resources import get_distribution, DistributionNotFound
#
# try:
#    # Change here if project is renamed and does not equal the package name
#    dist_name = "py-cordex"
#    __version__ = get_distribution(dist_name).version
# except DistributionNotFound:
#    __version__ = "unknown"
# finally:
#    del get_distribution, DistributionNotFound


from . import core, regions, tables, tutorial  # , cmor, preprocessing
from .core.domain import (
    cordex_domain,
    create_dataset,
    domain_info,
    rotated_coord_transform,
    vertices,
)
from .core.utils import map_crs
from .tables import domains, ecmwf
from .version import version

__version__ = version


__all__ = [
    "core",
    "regions",
    "tables",
    "tutorial",
    "cordex_domain",
    "create_dataset",
    "domain_info",
    "rotated_coord_transform",
    "vertices",
    "map_crs",
    "domains",
    "ecmwf",
]
