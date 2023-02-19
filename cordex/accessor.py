import xarray as xr

from .utils import _get_info, _guess_domain


class CordexAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    @property
    def domain_id(self):
        """Return the domain_id if possible."""
        return self._obj.attrs["CORDEX_domain"]

    @property
    def info(self):
        """Return domain info."""
        return _get_info(self._obj)

    @property
    def grid_mapping(self):
        return self._obj.cf["grid_mapping"]

    def guess(self):
        """Guess which domain."""
        return _guess_domain(self._obj)


#    def plot(self):
#        """Plot data on a map."""
#        return "plotting!"


@xr.register_dataset_accessor("cx")
class CordexDatasetAccessor(CordexAccessor):
    pass


@xr.register_dataarray_accessor("cx")
class CordexDataArrayAccessor(CordexAccessor):
    pass
