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


from . import core, regions#, cmor, preprocessing
from .core.domain import (
    cordex_domain,
    create_dataset,
    rotated_coord_transform,
    vertices,
    domain_info,
)

from .core.utils import map_crs

from .core import tutorial

from .tables import domains, ecmwf
from . import tables

from .version import version

__version__ = version
