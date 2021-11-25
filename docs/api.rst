.. currentmodule:: cordex

#############
API reference
#############

This page provides an auto-generated summary of the py-cordex API.


Top-level functions
===================

.. autosummary::
   :toctree: generated/

   cordex_domain
   domain_info
   create_dataset
   rotated_coord_transform
   map_crs
   vertices
   cordex_cmor_table

Tutorial
========

.. autosummary::
   :toctree: generated/

   tutorial.open_dataset
   tutorial.ensemble

Preprocessing
=============

.. autosummary::
   :toctree: generated/

   preprocessing.rename_cordex
   preprocessing.get_grid_mapping
   preprocessing.replace_coords
   preprocessing.cordex_dataset_id
   preprocessing.promote_empty_dims
   preprocessing.remap_lambert_conformal
   preprocessing.sort_ds_dict_by_attr
   preprocessing.get_grid_mapping_name
   preprocessing.member_id_to_dset_id
   preprocessing.attr_to_coord


Cmorization
===========

.. autosummary::
   :toctree: generated/

   cmor.season
   cmor.season_bounds
   cmor.mid_of_month
   cmor.month_bounds
   cmor.to_cftime

Regions within Euro-Cordex
==========================

.. autosummary::
   :toctree: generated/

   regions.germany
   regions.prudence
