from importlib.metadata import version as _get_version

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
from .tables import ecmwf
from .transform import (
    map_crs,
    rotated_coord_transform,
    transform,
    transform_bounds,
    transform_coords,
    derotate_vector,
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
    "ecmwf",
    "cell_area",
    "rewrite_coords",
    "domains",
]
