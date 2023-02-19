import xarray as xr

from .utils import _guess_domain


class CordexAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    @property
    def center(self):
        """Return the geographic center point of this dataset."""
        if self._center is None:
            # we can use a cache on our accessor objects, because accessors
            # themselves are cached on instances that access them.
            lon = self._obj.rlon
            lat = self._obj.rlat
            self._center = (float(lon.mean()), float(lat.mean()))
        return self._center

    @property
    def domain_id(self):
        """Return the geographic center point of this dataset."""
        return self._obj.attrs["CORDEX_domain"]

    @property
    def domain_info(self):
        """Return the geographic center point of this dataset."""
        return _guess_domain(self._obj)

    @property
    def grid_mapping(self):
        return self._obj.cf["grid_mapping"]

    def plot(self):
        """Plot data on a map."""
        return "plotting!"


@xr.register_dataset_accessor("cx")
class CordexDatasetAccessor(CordexAccessor):
    pass


@xr.register_dataarray_accessor("cx")
class CordexDataArrayAccessor(CordexAccessor):
    pass
