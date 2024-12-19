import pkg_resources

from . import regions, tables, tutorial
from .accessor import CordexDataArrayAccessor, CordexDatasetAccessor  # noqa
from .domain import (
    cordex_domain,
    create_dataset,
    domain,
    domain_info,
    vertices,
    rewrite_coords,
)
from .tables import domains, ecmwf
from .transform import (
    map_crs,
    rotated_coord_transform,
    transform,
    transform_bounds,
    transform_coords,
    derotate_vector,
)
from .utils import cell_area

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
    "domain",
    "cordex_domain",
    "create_dataset",
    "domain_info",
    "rotated_coord_transform",
    "vertices",
    "map_crs",
    "transform",
    "transform_coords",
    "transform_bounds",
    "derotate_vector",
    "domains",
    "ecmwf",
    "cell_area",
    "rewrite_coords",
]
