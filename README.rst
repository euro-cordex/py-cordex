=====================
Cordex Python Package
=====================

.. image:: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml/badge.svg
    :target: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml
    
.. image:: https://codecov.io/gh/euro-cordex/py-cordex/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/euro-cordex/py-cordex

.. image:: https://img.shields.io/pypi/v/py-cordex.svg
    :target: https://pypi.python.org/pypi/py-cordex
    
.. image:: https://anaconda.org/conda-forge/py-cordex/badges/installer/conda.svg
    :target: https://anaconda.org/conda-forge/py-cordex

.. image:: https://readthedocs.org/projects/py-cordex/badge/?version=latest
    :target: https://py-cordex.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://pyup.io/repos/github/euro-cordex/py-cordex/shield.svg
    :target: https://pyup.io/repos/github/euro-cordex/py-cordex/
    :alt: Updates

.. image:: https://zenodo.org/badge/304687410.svg
   :target: https://zenodo.org/badge/latestdoi/304687410

Python tools for the Cordex Community.

* Free software: MIT license
* Documentation: https://py-cordex.readthedocs.io.

Why
---

.. image:: http://mybinder.org/badge_logo.svg
    :alt: py-cordex examples
    :target: https://mybinder.org/v2/gh/euro-cordex/py-cordex/master?urlpath=lab%2Ftree%2Fnotebooks%2Fdomains.ipynb

Python tools to work with `CORDEX <https://cordex.org/>`_ domains, meta data and cmorization.

Meta Data
^^^^^^^^^
Access to meta data should be automatic and machine readable to avoid humans to do boring, repitive tasks that are error-prone. For that purpose, easy access to
meta information should be guaranteed by tables collected here: https://github.com/euro-cordex/tables 

Features
--------

* Tools to manage CORDEX grids as xarray datasets.
* Includes coordinate transformations, bounds and vertices for CORDEX datasets.
* Utitlities for cmorization to make the CORDEX ensembles more consistent.
* Preprocessing for easy access to a homogenized CORDEX ensemble dataset.

For planned features, please have a look at the `issues <https://github.com/euro-cordex/py-cordex/issues>`_, grab one, and `contribute <https://py-cordex.readthedocs.io/en/latest/contributing.html>`_!

Installation
------------

We recommend installing `py-cordex` with conda:

.. code-block:: console

    conda install -c conda-forge py-cordex
    

Installation from source
------------------------

You can install the package directly from github using pip:


.. code-block:: console

    pip install git+https://github.com/euro-cordex/py-cordex


If you want to contribute, please fork the repository to your github account
and install it in development mode, e.g.


.. code-block:: console

    git clone https://github.com/<your-account>/py-cordex
    cd cordex
    pip install -e .


This will install the package but you can still edit it and you don't need the package in your :code:`PYTHONPATH`.
Please get in contact as early as possible, e.g., using `draft pull requests <https://github.blog/2019-02-14-introducing-draft-pull-requests>`_.


Requirements
------------

* python3.6 or higher
* numpy
* pandas
* xarray
* netCDF4

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

Parts of this package have been developed within the project Pilot Lab Exascale Earth System Modelling (`PL-ExaESM <https://www.fz-juelich.de/SharedDocs/Meldungen/IAS/JSC/EN/2019/2019-09-pl-exaesm.html>`_).

