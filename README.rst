
py-cordex: create cordex grids and meta data
============================================

+----------------------------+-----------------------------------------------------+
| Versions                   | |pypi| |conda|                                      |
+----------------------------+-----------------------------------------------------+
| Documentation              | |docs| |versions| |binder|                          |
+----------------------------+-----------------------------------------------------+
| Open Source                | |license| |fair| |fossa| |zenodo|                   |
+----------------------------+-----------------------------------------------------+
| Coding Standards           | |black| |pre-commit|                                |
+----------------------------+-----------------------------------------------------+
| Development Status         | |ci| |codecov|                                      |
+----------------------------+-----------------------------------------------------+

This package offers python tools for the `CORDEX <https://cordex.org/>`_ community and should make your work with CORDEX grids and meta data easy.
Most of the tools leverage the ``xarray`` API to create grid and coordinate informations and data of CORDEX domains in the
form of an ``xarray.Dataset`` directly from the official `CORDEX archive specifications <https://cordex.org/experiment-guidelines/experiment-protocol-rcms/>`_.

Features
--------

* Tools to manage CORDEX grids as xarray datasets.
* Includes coordinate transformations, bounds and vertices for CORDEX datasets.
* Utitlities for cmorization to make the CORDEX ensembles more consistent.
* Preprocessing for easy access to a homogenized CORDEX ensemble dataset.

For planned features, please have a look at the `issues <https://github.com/euro-cordex/py-cordex/issues>`_, grab one, and `contribute <https://py-cordex.readthedocs.io/en/latest/contributing.html>`_!

Meta data
^^^^^^^^^
Access to meta data should be automatic and machine readable to avoid humans to do boring, repetitiv tasks that are error prone.
For that purpose, easy access to meta information should be guaranteed by tables collected `here <https://github.com/euro-cordex/tables>`_.

Installation
------------

We recommend installing `py-cordex` with conda:

.. code-block:: console

    conda install -c conda-forge py-cordex


Installation from source
------------------------

We don't recommend to pip install py-cordex because some of the dependencies require pre-compiled packages
that won't work with pip. For instructions to install py-cordex from source, please have a look
at the `contributing guide <https://py-cordex.readthedocs.io/en/stable/contributing.html>`_.
If you want to contribute, please get in contact as early as possible, e.g.,  using `draft pull requests <https://github.blog/2019-02-14-introducing-draft-pull-requests>`_.


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

Parts of this package have been developed within the project Pilot Lab Exascale Earth System Modelling (`PL-ExaESM <https://www.fz-juelich.de/SharedDocs/Meldungen/IAS/JSC/EN/2019/2019-09-pl-exaesm.html>`_).


.. |pypi| image:: https://img.shields.io/pypi/v/py-cordex.svg
        :target: https://pypi.python.org/pypi/py-cordex
        :alt: Python Package Index Build

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/py-cordex.svg
        :target: https://anaconda.org/conda-forge/py-cordex
        :alt: Conda-forge Build Version

.. |ci| image:: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml/badge.svg
        :target: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml
        :alt: Build Status

.. |codecov| image:: https://codecov.io/gh/euro-cordex/py-cordex/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/euro-cordex/py-cordex
        :alt: Covecov

.. |docs| image:: https://readthedocs.org/projects/py-cordex/badge
        :target: https://py-cordex.readthedocs.io/en/latest
        :alt: Documentation Status

.. |binder| image:: http://mybinder.org/badge_logo.svg
        :target: https://mybinder.org/v2/gh/WCRP-CORDEX/binder-sandbox/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252Feuro-cordex%252Fpy-cordex%26urlpath%3Dlab%252Ftree%252Fpy-cordex%252Fnotebooks%252Fdomains.ipynb%26branch%3Dmaster
        :alt: py-cordex examples

.. |zenodo| image:: https://zenodo.org/badge/304687410.svg
        :target: https://zenodo.org/badge/latestdoi/304687410
        :alt: DOI

.. |license| image:: https://img.shields.io/github/license/euro-cordex/py-cordex.svg
        :target: https://github.com/euro-cordex/py-cordex/blob/master/LICENSE
        :alt: License

.. |fair| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow
        :target: https://fair-software.eu
        :alt: FAIR Software Compliance

.. |fossa| image:: https://app.fossa.com/api/projects/git%2Bgithub.com%2Feuro-cordex%2Fpy-cordex.svg?type=shield
        :target: https://app.fossa.com/projects/git%2Bgithub.com%2Feuro-cordex%2Fpy-cordex?ref=badge_shield
        :alt: FOSSA

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black
        :alt: Python Black

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/euro-cordex/py-cordex/master.svg
        :target: https://results.pre-commit.ci/latest/github/euro-cordex/py-cordex/master
        :alt: pre-commit.ci status

.. |versions| image:: https://img.shields.io/pypi/pyversions/py-cordex.svg
        :target: https://pypi.python.org/pypi/py-cordex
        :alt: Supported Python Versions

.. |funding| image:: https://img.shields.io/badge/Powered%20by-ExaESM-blue.svg
        :target: https://www.exaesm.de/
        :alt: Funding
