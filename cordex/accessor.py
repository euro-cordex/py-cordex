import xarray as xr

from .utils import _get_info, _guess_domain


class CordexAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    @property
    def domain_id(self):
        """Returns the domain_id."""
        return self._obj.attrs["CORDEX_domain"]

    @property
    def grid_mapping(self):
        """Returns the grid_mapping variable."""
        return self._obj.cf["grid_mapping"]

    def info(self):
        """Return domain info in CORDEX format.

        The function returns a dictionary containing domain
        information in the format of the CORDEX archive specifications.

        Returns
        -------
        domain info : dict
            A dictionary that contains domain information.

        """
        return _get_info(self._obj)

    def guess(self):
        """Guess which domain this could be.

        Compares the coordinate axis information to known the
        coordinates of known CORDEX domains to guess the
        ``domain_id``.

        Returns
        -------
        domain info : dict
            A dictionary that contains domain information.

        """
        return _guess_domain(self._obj)


@xr.register_dataset_accessor("cx")
class CordexDatasetAccessor(CordexAccessor):
    pass


@xr.register_dataarray_accessor("cx")
class CordexDataArrayAccessor(CordexAccessor):
    pass
