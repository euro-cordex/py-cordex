# -*- coding: utf-8 -*-

import pkg_resources

from . import core, regions, tables, tutorial  # , cmor, preprocessing
from .core.domain import (
    cordex_domain,
    create_dataset,
    domain_info,
    rotated_coord_transform,
    vertices,
)
from .core.transform import map_crs, transform
from .tables import domains, ecmwf

# from .version import version

# __version__ = version

try:
    __version__ = pkg_resources.get_distribution("py-cordex").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"


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
    "transform",
    "domains",
    "ecmwf",
]
