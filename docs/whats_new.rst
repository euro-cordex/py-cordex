.. currentmodule:: cordex

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo

v0.4.0 (29 April 2022)
----------------------

This version introduces CORDEX regular grids, they are denoted with an i in the end, e.g., ``EUR-11i``. The cmor tables have also been
updated to allow for a preliminary cmorization of CORDEX data using CMIP6 cmor tables. However, this is still only used for
development since there is no official data request yet for downscaling CMIP6 models. There are also some internal changes,
mainly to improve code quality for easier contributions.

New Features
~~~~~~~~~~~~

- Updated ``map_crs`` to work for arbitray transformations (:pull:`52`).
- Added ``cordex-regular`` table and options to create regular cordex domains (:pull:`53`, :pull:`54`).

Internal Changes
~~~~~~~~~~~~~~~~

- Updated documentation including contribution guide, switch to ``sphinx-book-theme`` (:pull:`57`).
- Introduced ``linting.yaml`` test and ``pre-commit`` hook (:pull:`56`).
- Updated external resource from pooch to ignore hashes and use main branch from `tables <https://github.com/euro-cordex/tables>`_ repository (:pull:`53`). This should protect earlier version from breaking if the tables change.
- Changed ``pooch`` resource for CMOR table fetching. The test workflow now uses original `CMIP6 cmor tables <https://github.com/PCMDI/cmip6-cmor-tables>`_ in a combination with an updated controlled vocabulary for `CORDEX_CV.json <https://github.com/euro-cordex/cordex-cmor-tables>`_. Because there are no official CORDEX CMIP6 CMOR tables yet, we ignore hash checking for now. This will change in the future (:pull:`55`).

Breaking Changes
~~~~~~~~~~~~~~~~

- ``map_crs`` switched order of coordinates to `COARDS <https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions>`_ (:pull:`52`).


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
