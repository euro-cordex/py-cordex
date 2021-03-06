=====================
Cordex Python Package
=====================


.. image:: https://img.shields.io/pypi/v/py-cordex.svg
    :target: https://pypi.python.org/pypi/py-cordex
        
.. image:: https://travis-ci.com/euro-cordex/py-cordex.svg?branch=master
    :target: https://travis-ci.com/euro-cordex/py-cordex

.. image:: https://readthedocs.org/projects/py-cordex/badge/?version=latest
    :target: https://py-cordex.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/euro-cordex/py-cordex/badge.svg?branch=master
    :target: https://coveralls.io/github/euro-cordex/py-cordex?branch=master

.. image:: https://pyup.io/repos/github/euro-cordex/py-cordex/shield.svg
    :target: https://pyup.io/repos/github/euro-cordex/py-cordex/
    :alt: Updates



Python tools for the Cordex Community.


* Free software: MIT license
* Documentation: https://py-cordex.readthedocs.io.

Why
---

The goal of this software package is to have a unified common code base for typical and reacurring tasks in the `CORDEX <https://cordex.org/>`_ community.
More precisely, there are several tasks related to regional climate model data processing that are model independent and could become unified. This should
enable better reproducability and exchange of scientific analyses in the community.

Meta Data
^^^^^^^^^
Access to meta data should be automatic and machine readable to avoid humans to do boring, repitive tasks that are error-prone. For that purpose, easy access to
meta information should be guaranteed by tables collected here: https://github.com/euro-cordex/tables 

Features
--------

* Tools to manage Cordex Domain Data and Datasets
* Easy Cordex ESGF Access

For planned features, please have a look at the `issues <https://github.com/euro-cordex/py-cordex/issues>`_, grab one, and `contribute <https://py-cordex.readthedocs.io/en/latest/contributing.html>`_!

Installation
------------

You can install the package directly from github using pip:


.. code-block:: console

    pip install git+https://github.com/euro-cordex/py-cordex


If you want to contribute, I recommend cloning the repository and installing the package in development mode, e.g.


.. code-block:: console

    git clone https://github.com/euro-cordex/py-cordex
    cd cordex
    pip install -e .


This will install the package but you can still edit it and you don't need the package in your :code:`PYTHONPATH`

Highlights
----------

Use the domain module to create Cordex domains safe and easy:

.. code-block:: python

    from cordex import domain as dm
    eur11 = domain('EUR-11')
    

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
