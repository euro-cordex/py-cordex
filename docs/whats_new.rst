.. currentmodule:: cordex

What's New
==========

.. ipython:: python
   :suppress:

    import cordex

v0.5.0 (2 March 2023)
---------------------

This release introduces a cmorization workflow using :py:meth:`cmor.cmorize_variable` and deprecates some code depending on `cartopy`.
There has also been a lot of refactoring since `py-cordex` now makes extensive use of `cf-xarray <https://github.com/xarray-contrib/cf-xarray>`_
to make cmorization as easy as possible.
This release is used for developing and testing the new `CORDEX cmor tables <https://github.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/tree/main/Tables>`_
for downscaling CMIP6 global model runs. Since the new archive specifications for CORDEX are not fixed yet, the API will probably change in the future.
For some examples of cmorization, please visit the `CORDEX cmorization examples notebook <https://wcrp-cordex.github.io/cordex-cmip6-cmor-tables/cmor-examples.html>`_ in the cmor table repository. This version also drops `python3.7` support.

New Features
~~~~~~~~~~~~

- Cmorization module from ``pyremo.cmor`` has been moved upstream into ``cordex.cmor``. A new cmorization funcion :py:meth:`cmor.cmorize_variable` now allows for easy cmorization using the archive specifications included in ``py-cordex`` (:pull:`68`, :pull:`74`, :pull:`76`, :pull:`77`, :pull:`78`, :pull:`92`).

- Added :py:meth:`transform`, :py:meth:`transform_coords` and :py:meth:`transform_bounds` in favour of deprecated :py:meth:`map_crs`, :py:meth:`rotated_coord_transform` and :py:meth:`vertices`. Dropped ``cartopy`` dependency in favour of ``pyproj`` (:pull:`71`, :pull:`102`).

- Introduced ``cx`` Dataset and DataArray accessors (:pull:`105`, :pull:`112`).

Breaking Changes
~~~~~~~~~~~~~~~~

* The standard cmor tables are now fetched from the `development repository <https://github.com/WCRP-CORDEX/cordex-cmip6-cmor-tables/tree/main/Tables>`_ (:pull:`116`).

* ``nan`` values in domain tables are replace with ``None``. That means in case you retrive domain information with :py:meth:`domain_info` for a regular grid, ``pollon`` and ``pollat`` will now be ``None`` instead of ``nan`` (:issue:`107`, :pull:`105`).

Deprecations
~~~~~~~~~~~~

- ``add_vertices`` keyword is deprecated in favour of ``bounds`` keyword in :py:meth:`cordex_domain` and :py:meth:`create_dataset`. Note that the 2D longitude and latitude bounds (``lon_vertices``, ``lat_vertices``) will now appear as coordinates in the domain dataset by default and not as data variables (:pull:`101`).
- :py:meth:`map_crs` and :py:meth:`rotated_coord_transform` are deprecated and will be removed in the future (:pull:`71`).
- :py:meth:`vertices` is deprecated, please use :py:meth:`transform_bounds` instead (:pull:`102`).
- :py:meth:`tables.cmip6_cmor_table` is deprecated in favour of :py:meth:`tables.cordex_cmor_table` (:pull:`117`).

Internal Changes
~~~~~~~~~~~~~~~~

- Updated CI pipeline (:pull:`83`).
- Refactored transformation functions (:pull:`71`).
- Pinned ``pyproj>=3.3.0`` and ``cf_xarray>=0.8.0`` (Dropped python 3.7 support).
- Moved core modules (:pull:`87`).
- Added ``CITATION.cff`` (:pull:`87`).
- Updated documentation environment for python 3.10 (:pull:`89`).
- Use ``cf-xarray`` for adding bounds to regular datasets (:pull:`101`).



v0.4.1 (23 June 2022)
---------------------

Patch release to fix coordinate precisions, also includes some maintenance updates.

Internal Changes
~~~~~~~~~~~~~~~~

- Added pre-commit hooks for Jupyter notebooks (including reformatting), added ``.zenodo.json`` (:pull:`59`).
- Switched to automatic version numbering using ``setuptools_scm`` (:pull:`62`).
- Added ``publish-pypi.yaml`` workflow (:pull:`63`).

Bug Fixes
~~~~~~~~~

- Fixed rounding errors for grid coordinates (:pull:`61`). This might have some influence on results, since coordinate values might change to more exact values.


v0.4.0 (29 April 2022)
----------------------

This version introduces CORDEX regular grids, they are denoted with an i in the end, e.g., ``EUR-11i``. The cmor tables have also been
updated to allow for a preliminary cmorization of CORDEX data using CMIP6 cmor tables. However, this is still only used for
development since there is no official data request yet for downscaling CMIP6 models. There are also some internal changes,
mainly to improve code quality for easier contributions.

New Features
~~~~~~~~~~~~

- Updated :py:meth:`map_crs` to work for arbitray transformations (:pull:`52`).
- Added ``cordex-regular`` table and options to create regular cordex domains (:pull:`53`, :pull:`54`).

Internal Changes
~~~~~~~~~~~~~~~~

- Updated documentation including contribution guide, switch to ``sphinx-book-theme`` (:pull:`57`).
- Introduced ``linting.yaml`` test and ``pre-commit`` hook (:pull:`56`).
- Updated external resource from pooch to ignore hashes and use main branch from `tables <https://github.com/euro-cordex/tables>`_ repository (:pull:`53`). This should protect earlier version from breaking if the tables change.
- Changed ``pooch`` resource for CMOR table fetching. The test workflow now uses original `CMIP6 cmor tables <https://github.com/PCMDI/cmip6-cmor-tables>`_ in a combination with an updated controlled vocabulary for `CORDEX_CV.json <https://github.com/euro-cordex/cordex-cmor-tables>`_. Because there are no official CORDEX CMIP6 CMOR tables yet, we ignore hash checking for now. This will change in the future (:pull:`55`).

Breaking Changes
~~~~~~~~~~~~~~~~

- :py:meth:`map_crs` switched order of coordinates to `COARDS <https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions>`_ (:pull:`52`).


v0.3.2 (30 March 2022)
----------------------

Patch release to fix vertices issues, imports and dependencies. Also includes ``python3.10`` support.

Internal Changes
~~~~~~~~~~~~~~~~

- Added ``python3.10`` support (:pull:`48`).
- Fixed vertices attributes (:pull:`49`).

Breaking Changes
~~~~~~~~~~~~~~~~

- Fixed vertices order (:pull:`50`).
- Fixed imports (:pull:`51`).


v0.3.1 (24 February 2022)
-------------------------

Patch release to fix dependencies when ``py-cordex`` is installed with pip.

Bug Fixes
~~~~~~~~~
Added ``cftime`` to pip dependencies.


v0.3.0 (12 January 2022)
------------------------

This release introduces the cmorization and preprocessing modules. We outsourced
some more code that was REMO dependent to be more useful for CORDEX datasets
in general. The CORDEX tables are not available for CMIP6 yet, so the tables
are only useful for the development of a cmorization workflow. The peprocessing
module aims at overcoming CF related inconsistencies in ESGF CORDEX datasets.

New Features
~~~~~~~~~~~~
- Included experiment cmor table resource for development (:pull:`27`).
- Included basic utitilities to compute time bounds (:pull:`28`).
- Added ``add_vertices`` option to ``cx.cordex_domain`` (:pull:`32`).
- Additional cmorization utilities (:pull:`28`, :pull:`30`).
- Preprocessing module for working with CORDEX ensembles (:pull:`35`, :pull:`36`).
- Tutorial module for loading simple CORDEX test datasets (:pull:`35`).

Internal Changes
~~~~~~~~~~~~~~~~
- Coverage is ignored for unused modules (:pull:`29`).


v0.2.1 (1 November 2021)
------------------------

Small bugfix release that updates the path to download germany shapefiles (:pull:`24`).
Also fixed coordinate attribute issues (:pull:`25`) and added a test.


v0.2.0 (28 October 2021)
------------------------

This is a new restructuring release, that cleanded a lot of things and includes
much more documentation and example notebooks. All meta information is removed
and accessed online. This version should be more open for contributions!

New Features
~~~~~~~~~~~~
- Included new sub regions (germany and prudence) for masking and analysis.
- Included function ``map_crs`` for coordinate transformations using cartopy.

Internal Changes
~~~~~~~~~~~~~~~~
- Tables are removed from the package and stored in an extra github repo.
- Tables are download at first access using pooch.
- Using now conda environemt for testing (:pull:`18`).
- Restructured and cleaned dependencies (:pull:`21`).


v0.1.2 (3 June 2021)
--------------------

This is a major restructuring release. The code base has been reduced significantly
and the main data structure are now xarray's datarrays. Several new domains have been
added.
