import pkg_resources

from . import core, regions, tables, tutorial
from .core.domain import cordex_domain, create_dataset, domain_info, vertices
from .core.transform import (
    map_crs,
    rotated_coord_transform,
    transform,
    transform_coords,
)
from .tables import domains, ecmwf

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
    "transform_coords",
    "domains",
    "ecmwf",
]
