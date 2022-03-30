.. currentmodule:: cordex

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo
    
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
