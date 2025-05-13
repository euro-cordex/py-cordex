import numpy as np
import xarray as xr

from .config import nround

# from .utils import _get_info, _guess_domain
from .tables import domains
from .domain import rewrite_coords

IDS = ["domain_id", "CORDEX_domain"]


def _get_domain_id(ds):
    """search for any valid domain id"""
    for attr in IDS:
        if attr in ds.attrs:
            return ds.attrs[attr]
    return None


def _get_info(ds, tables=None, precision=nround):
    try:
        grid_mapping_name = ds.cf["grid_mapping"].grid_mapping_name
    except KeyError:
        grid_mapping_name = None

    if grid_mapping_name and grid_mapping_name != "rotated_latitude_longitude":
        raise ValueError(
            f"Grid mapping name '{grid_mapping_name}' is not supported. Only 'rotated_latitude_longitude' is supported."
        )
    # check if tables is None
    if tables is None:
        tables = domains.table.replace(np.nan, None)
    try:
        x = ds.cf["X"]
        y = ds.cf["Y"]
    except KeyError:
        x = ds.rlon
        y = ds.rlat
    nlon = x.size
    nlat = y.size
    dlon = x.diff(x.dims[0]).values.item(0)
    dlat = y.diff(y.dims[0]).values.item(0)
    ll_lon = x[0].min().values.item(0)
    ll_lat = y[0].min().values.item(0)
    ur_lon = x[-1].min().values.item(0)
    ur_lat = y[-1].min().values.item(0)
    try:
        pollon = ds.cf["grid_mapping"].grid_north_pole_longitude
        pollat = ds.cf["grid_mapping"].grid_north_pole_latitude
    except KeyError:
        pollon = None
        pollat = None
    coords = {
        "nlon": nlon,
        "nlat": nlat,
        "ll_lon": ll_lon,
        "ur_lon": ur_lon,
        "ll_lat": ll_lat,
        "ur_lat": ur_lat,
        "dlon": dlon,
        "dlat": dlat,
        "pollon": pollon,
        "pollat": pollat,
    }
    # round
    info = {
        k: (np.round(v, nround) if isinstance(v, float) else v)
        for k, v in coords.items()
    }
    return info


def _guess_domain(ds, tables=None):
    if tables is None:
        tables = domains.table  # .replace(np.nan, None)
    try:
        info = _get_info(ds, tables)
    except Exception as e:
        print(e)
        raise Exception(
            "Could not determine domain, only rotated_latitude_longitude supported."
        )
    filt = tables
    for k, v in info.items():
        if filt.empty:
            return None
        if v:
            filt = filt[np.isclose(filt[k], v)]
        else:
            filt = filt[np.isnan(filt[k])]  # | filt[k] is None]
    # reset index and convert to dict
    return filt.reset_index().iloc[0].replace(np.nan, None).to_dict()


class CordexAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._domain_id = None
        self._info = None
        self._guess = None

    @property
    def domain_id(self, guess=True):
        """
        Returns the domain_id.

        This property will return the ``CORDEX_domain`` or ``domain_id`` global
        attribute if present. If none of those attributes are found, the
        domain information will be guessed.

        Parameters
        ----------
        guess : bool, optional
            If True, the domain information will be guessed if not found. Default is True.

        Returns
        -------
        str
            The domain_id.
        """
        if self._domain_id is None:
            self._domain_id = _get_domain_id(self._obj)
        if self._domain_id is None:
            self._domain_id = self.guess()["short_name"]
        return self._domain_id

    @property
    def grid_mapping(self):
        """
        Returns the grid_mapping variable.

        Returns
        -------
        xarray.DataArray
            The grid_mapping variable from the xarray object.
        """
        return self._obj.cf["grid_mapping"]

    def info(self):
        """
        Return domain info in CORDEX format.

        The function returns a dictionary containing domain
        information in the format of the CORDEX archive specifications.

        Returns
        -------
        dict
            A dictionary that contains domain information.
        """
        if self._info is None:
            self._info = _get_info(self._obj)
        return self._info

    def guess(self):
        """
        Guess which domain this could be.

        Compares the coordinate axis information to known
        coordinates of known CORDEX domains to guess the
        ``domain_id``.

        Returns
        -------
        dict
            A dictionary containing the guessed domain information.
        """
        if self._guess is None:
            self._guess = _guess_domain(self._obj)
        return self._guess

    def map(self, projection=None):
        """Create a simple overview map.

        Creates a simple map overview. By default, the map projection
        defaults to the grid mapping attribute.

        Parameters
        ----------
        projection : cartopy.crs
            CRS used for projection. By default, the map projection
            defaults to the CRS defined by the grid mapping attribute.

        Returns
        -------
        ax : GeoAxesSubplot
            Cartopy plot projection using tiles.

        """
        import cartopy.crs as ccrs
        import cartopy.feature as cf

        # import cartopy.io.img_tiles as cimgt
        import matplotlib.pyplot as plt

        obj = self._obj

        try:
            x = obj.cf["X"]
            y = obj.cf["Y"]
        except KeyError:
            x = obj.rlon
            y = obj.rlat

        mapping = obj.cf["grid_mapping"]
        pole = (
            mapping.grid_north_pole_longitude,
            mapping.grid_north_pole_latitude,
        )
        central_longitude = 0.0
        if x.min() > 0.0:
            central_longitude = 180.0

        transform = ccrs.RotatedPole(*pole, central_rotated_longitude=central_longitude)

        if projection is None:
            projection = transform

        ax = plt.axes(projection=projection)
        # use google maps tiles
        ax.stock_img()
        # request = cimgt.GoogleTiles()
        # ax.add_image(request, 4)  # , interpolation='spline36', regrid_shape=2000)
        ax.gridlines(
            draw_labels=True,
            linewidth=0.5,
            color="gray",
            xlocs=range(-180, 180, 10),
            ylocs=range(-90, 90, 5),
        )
        ax.set_extent(
            [
                x.min() - central_longitude,
                x.max() - central_longitude,
                y.min(),
                y.max(),
            ],
            crs=transform,
        )
        ax.coastlines(resolution="110m", color="black", linewidth=1)
        ax.add_feature(cf.BORDERS, color="black")

        return ax
        # ax.set_title(CORDEX_domain)

    def rewrite_coords(
        self,
        coords="xy",
        bounds=False,
        domain_id=None,
        mip_era="CMIP5",
        method="nearest",
    ):
        """
        Rewrite coordinates in a dataset.

        This function ensures that the coordinates in a dataset are consistent and can be
        compared to other datasets. It can reindex the dataset based on specified coordinates
        or domain information while trying to keep the original coordinate attributes.

        Parameters
        ----------
        coords : str, optional
            Specifies which coordinates to rewrite. Options are:
            - "xy": Rewrite only the X and Y coordinates.
            - "lonlat": Rewrite only the longitude and latitude coordinates.
            - "all": Rewrite both X, Y, longitude, and latitude coordinates.
            Default is "xy". If longitude and latitude coordinates are not present in the dataset, they will be added.
            Rewriting longitude and latitude coordinates is only possible if the dataset contains a grid mapping variable.
        bounds : bool, optional
            If True, the function will also handle the bounds of the coordinates. If the dataset already has bounds,
            they will be updated while preserving attributes and shape. If not, the bounds will be assigned.
        domain_id : str, optional
            The domain identifier used to obtain grid information. If not provided, the function will attempt
            to use the domain_id attribute from the dataset.
        mip_era : str, optional
            The MIP era (e.g., "CMIP5", "CMIP6") used to determine coordinate attributes. Default is "CMIP5".
            Only used if the dataset does not already contain coordinate attributes.
        method : str, optional
            The method used for reindexing the X and Y axis. Options include "nearest", "linear", etc. Default is "nearest".

        Returns
        -------
        ds : xr.Dataset
            The dataset with rewritten coordinates.
        """
        return rewrite_coords(
            self._obj,
            coords=coords,
            bounds=bounds,
            domain_id=domain_id,
            mip_era=mip_era,
            method=method,
        )


@xr.register_dataset_accessor("cx")
class CordexDatasetAccessor(CordexAccessor):
    pass


@xr.register_dataarray_accessor("cx")
class CordexDataArrayAccessor(CordexAccessor):
    pass
