=====================
Cordex Python Package
=====================

.. image:: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml/badge.svg
    :target: https://github.com/euro-cordex/py-cordex/actions/workflows/ci.yaml
    
.. image:: https://codecov.io/gh/euro-cordex/py-cordex/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/euro-cordex/py-cordex

.. image:: https://img.shields.io/pypi/v/py-cordex.svg
    :target: https://pypi.python.org/pypi/py-cordex

.. image:: https://readthedocs.org/projects/py-cordex/badge/?version=latest
    :target: https://py-cordex.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://pyup.io/repos/github/euro-cordex/py-cordex/shield.svg
    :target: https://pyup.io/repos/github/euro-cordex/py-cordex/
    :alt: Updates



Python tools for the Cordex Community.

* Free software: MIT license
* Documentation: https://py-cordex.readthedocs.io.

Why
---

Python tools to work with `CORDEX <https://cordex.org/>`_ domains and meta data.

Meta Data
^^^^^^^^^
Access to meta data should be automatic and machine readable to avoid humans to do boring, repitive tasks that are error-prone. For that purpose, easy access to
meta information should be guaranteed by tables collected here: https://github.com/euro-cordex/tables 

Features
--------

* Tools to manage Cordex Domain Data and Datasets

For planned features, please have a look at the `issues <https://github.com/euro-cordex/py-cordex/issues>`_, grab one, and `contribute <https://py-cordex.readthedocs.io/en/latest/contributing.html>`_!

Installation from source
------------------------

You can install the package directly from github using pip:


.. code-block:: console

    pip install git+https://github.com/euro-cordex/py-cordex


If you want to contribute, I recommend cloning the repository and installing the package in development mode, e.g.


.. code-block:: console

    git clone https://github.com/euro-cordex/py-cordex
    cd cordex
    pip install -e .


This will install the package but you can still edit it and you don't need the package in your :code:`PYTHONPATH`


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
