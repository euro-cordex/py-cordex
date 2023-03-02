import xarray as xr

from .utils import _get_info, _guess_domain

IDS = ["domain_id", "CORDEX_domain"]


def _get_domain_id(ds):
    """search for any valid domain id"""
    for attr in IDS:
        if attr in ds.attrs:
            return ds.attrs[attr]
    return None


class CordexAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._domain_id = None
        self._info = None
        self._guess = None

    @property
    def domain_id(self, guess=True):
        """Returns the domain_id.

        This property will return the ``CORDEX_domain`` or ``domain_id` global
        attribute if present. If none of those attributes are found, the
        domain information will be guessed.

        """
        if self._domain_id is None:
            self._domain_id = _get_domain_id(self._obj)
        if self._domain_id is None:
            self._domain_id = self.guess()["short_name"]
        return self._domain_id

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
        if self._info is None:
            self._info = _get_info(self._obj)
        return self._info

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
        if self._guess is None:
            self._guess = _guess_domain(self._obj)
        return self._guess


@xr.register_dataset_accessor("cx")
class CordexDatasetAccessor(CordexAccessor):
    pass


@xr.register_dataarray_accessor("cx")
class CordexDataArrayAccessor(CordexAccessor):
    pass
