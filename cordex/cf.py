# -*- coding: utf-8 -*-
# flake8: noqa

import numpy as np

LAT_NAME = "lat"
LON_NAME = "lon"
RLAT_NAME = "rlat"
RLON_NAME = "rlon"
LON_BOUNDS = "lon_vertices"
LAT_BOUNDS = "lat_vertices"
BOUNDS_DIM = "vertices"

coords = {
    "rlon": {
        "axis": "X",
        "standard_name": "grid_longitude",
        "long_name": "longitude in rotated pole grid",
        "units": "degrees",
    },
    "rlat": {
        "axis": "Y",
        "standard_name": "grid_latitude",
        "long_name": "latitude in rotated pole grid",
        "units": "degrees",
    },
    "lon": {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
    },
    "lat": {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
    },
    "lon_vertices": {
        "units": "degrees_east",
    },
    "lat_vertices": {
        "units": "degrees_north",
    },
}


mapping = {
    "grid_mapping_name": "rotated_latitude_longitude",
    "grid_north_pole_latitude": 90.0,
    "grid_north_pole_longitude": 180.0,
    "north_pole_grid_longitude": 0.0,
}


DEFAULT_MAPPING_NCVAR = "rotated_latitude_longitude"


DEFAULT_CORDEX_ATTRS = {
    "institution": "",
    "institute_id": "",
    "experiment_id": "",
    "source": "",
    "model_id": "",
    "contact": "",
    "comment": "",
    "references": "",
    "initialization_method": 0,
    "physics_version": 0,
    "tracking_id": "",
    "CORDEX_domain": "",
    "driving_experiment": ", , ",
    "driving_model_id": "",
    "driving_model_ensemble_member": "r0i0p0",
    "driving_experiment_name": "",
    "rcm_version_id": "v1",
    "product": "output",
    "experiment": "",
    "frequency": "",
    "creation_date": "",
    "Conventions": "CF-1.4",
    "project_id": "CORDEX",
    "table_id": "",
    "title": "",
    "modeling_realm": "atmos",
    "realization": 0,
    "cmor_version": "",
}

grid_mapping_dtype = np.int32


vocabulary = {
    "CMIP5": dict(
        domain_id="CORDEX_domain",
        dims={
            "LAT": "lat",
            "LON": "lon",
            "Y": "rlat",
            "X": "rlon",
            "LON_BOUNDS": "lon_vertices",
            "LAT_BOUNDS": "lat_vertices",
            "BOUNDS_DIM": "vertices",
        },
        coords={
            "rlon": {
                "axis": "X",
                "standard_name": "grid_longitude",
                "long_name": "longitude in rotated pole grid",
                "units": "degrees",
            },
            "rlat": {
                "axis": "Y",
                "standard_name": "grid_latitude",
                "long_name": "latitude in rotated pole grid",
                "units": "degrees",
            },
            "lon": {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
            },
            "lat": {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
            },
            "lon_vertices": {
                "units": "degrees_east",
            },
            "lat_vertices": {
                "units": "degrees_north",
            },
        },
        default_mapping_ncvar="rotated_latitude_longitude",
        default_global_attrs={
            "institution": "",
            "institute_id": "",
            "experiment_id": "",
            "source": "",
            "model_id": "",
            "contact": "",
            "comment": "",
            "references": "",
            "initialization_method": 0,
            "physics_version": 0,
            "tracking_id": "",
            "CORDEX_domain": "",
            "driving_experiment": ", , ",
            "driving_model_id": "",
            "driving_model_ensemble_member": "r0i0p0",
            "driving_experiment_name": "",
            "rcm_version_id": "v1",
            "product": "output",
            "experiment": "",
            "frequency": "",
            "creation_date": "",
            "Conventions": "CF-1.4",
            "project_id": "CORDEX",
            "table_id": "",
            "title": "",
            "modeling_realm": "atmos",
            "realization": 0,
            "cmor_version": "",
        },
    ),
    "CMIP6": dict(
        domain_id="domain_id",
        dims={
            "LAT": "lat",
            "LON": "lon",
            "Y": "rlat",
            "X": "rlon",
            "LON_BOUNDS": "vertices_lon",
            "LAT_BOUNDS": "vertices_lat",
            "BOUNDS_DIM": "vertices",
        },
        coords={
            "rlon": {
                "axis": "X",
                "standard_name": "grid_longitude",
                "long_name": "longitude in rotated pole grid",
                "units": "degrees",
            },
            "rlat": {
                "axis": "Y",
                "standard_name": "grid_latitude",
                "long_name": "latitude in rotated pole grid",
                "units": "degrees",
            },
            "lon": {
                "standard_name": "longitude",
                "long_name": "longitude",
                "units": "degrees_east",
            },
            "lat": {
                "standard_name": "latitude",
                "long_name": "latitude",
                "units": "degrees_north",
            },
            "vertices_lon": {
                "units": "degrees_east",
            },
            "vertices_lat": {
                "units": "degrees_north",
            },
        },
        default_mapping_ncvar="crs",
        default_global_attrs={
            "institution": "",
            "institution_id": "",
            "experiment_id": "",
            "source_id": "",
            "contact": "",
            "comment": "",
        },
    ),
}
