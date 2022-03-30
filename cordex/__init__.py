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
    map_crs,
    vertices,
    domain_info,
)

from .core import tutorial

from .tables import domains, ecmwf, cordex_cmor_table
from . import tables

from .version import version

__version__ = version
__cmor_table_version__ = tables.__cmor_table_version__
