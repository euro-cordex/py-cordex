from importlib.metadata import version as _get_version

from . import regions, tables, tutorial
from .accessor import CordexDataArrayAccessor, CordexDatasetAccessor  # noqa
from .domain import (
    cordex_domain,
    create_dataset,
    domain,
    domain_info,
    rewrite_coords,
    vertices,
)
from .tables import ecmwf
from .transform import (
    derotate_vector,
    map_crs,
    rotated_coord_transform,
    transform,
    transform_bounds,
    transform_coords,
)
from .utils import cell_area


# keep this for backward compatibility
class domains:
    table = tables.domains


try:
    __version__ = _get_version("py-cordex")
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"


__all__ = [
    "cell_area",
    "cordex_domain",
    "core",
    "create_dataset",
    "derotate_vector",
    "domain",
    "domain_info",
    "domains",
    "ecmwf",
    "map_crs",
    "regions",
    "rewrite_coords",
    "rotated_coord_transform",
    "tables",
    "transform",
    "transform_bounds",
    "transform_coords",
    "tutorial",
    "vertices",
]
