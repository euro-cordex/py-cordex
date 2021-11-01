.. currentmodule:: cordex

What's New
==========

.. ipython:: python
   :suppress:

    import pyremo 

v0.2.1 (01 November 2021)
-------------------------

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
