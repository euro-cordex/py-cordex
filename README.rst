=====================
Cordex Python Package
=====================


.. image:: https://img.shields.io/pypi/v/py-cordex.svg
        :target: https://pypi.python.org/pypi/py-cordex

.. image:: https://img.shields.io/travis/euro-cordex/cordex.svg
        :target: https://travis-ci.com/euro-cordex/cordex

.. image:: https://readthedocs.org/projects/cordex/badge/?version=latest
        :target: https://cordex.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/euro-cordex/cordex/shield.svg
     :target: https://pyup.io/repos/github/euro-cordex/cordex/
     :alt: Updates



Python tools for the Cordex Community.


* Free software: MIT license
* Documentation: https://cordex.readthedocs.io.


Features
--------

* Manage Cordex Meta Data
* Tools to manage Cordex Domain Data
* ...

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

    # available tables
    print(dm.tables())
    # available domains names
    print(dm.names())
    # available domains names
    print(dm.names('cordex-core'))
    # print cordex core table
    print(dm.table('cordex-core'))

    # show domains with some dummy data (uses cdo topo)
    for short_name, domain in dm.domains().items():
        print(domain)
        domain.to_netcdf(short_name+'.nc', dummy='topo')

Using the ESGF module for browsing the file system and creating pandas dataframes of file selections:
(DKRZ example).

    Note that DKRZ also provides an **intake-esm** table (https://intake-esm.readthedocs.io/en/latest/) 
    for browsing and searchin MIP data at DKRZ. I highly recommend to have a look at it and create tables
    for your own data pool. For creating a table, you could use our cordex conventions module!

.. code-block:: python

    import datetime as dt
    from cordex.ESGF import get_selection, conventions

    project_id = 'CMIP5'
    root       = '/pool/data/CMIP5/cmip5'

    filter   = {'institute'       : 'MPI-M',
                'product'         : 'output1',
                'variable'        : 'tas',
                'model'           : 'MPI-ESM-LR',
                'modeling_realm'  : 'atmos',
                'experiment'      : 'historical',
                'ensemble_member' : 'r1i1p1'}

    # get a selection of files (using pandas dataframes)
    selection = get_selection(project_id, root=root, filter=filter)
    print(selection)

    # create a finer selection and convert dates to datetime objects
    selection = selection.subset(variable='pr').to_datetime()
    # get a timeseries of files
    selection = selection.select_timerange([dt.datetime(1990,1,1),dt.datetime(2000,1,1)])
    print(selection)

Use the ESGF module to create your filename using an attribute dictionary. Use the CORDEX filenaming
convetion (or create your own!):

.. code-block:: python

    from cordex import ESGF as esgf

    root       = '/my_root'

    # define attributes
    attributes   = {'institute_id'    : 'GERICS',
                    'product'         : 'output',
                    'model_id'        : 'GERICS-REMO2015',
                    'experiment_id'   : 'evaluation',
                    'driving_model_id': 'ECMWF-ERAINT',
                    'variable'        : 'pr',
                    'rcm_version_id'  : 'v1',
                    'date'            : 'v20200221',
                    'frequency'       : 'day',
                    'CORDEX_domain'   : 'EUR-11',
                    'suffix'          : 'nc',
                    'ensemble'        : 'r1i1p1'}

    # we use the CORDEX convention as example
    convention = esgf.CORDEX()
    # print the convention patterns 
    print(convention.path_conv.conv_str)
    print(convention.filename_conv.conv_str)
    # only filename
    filename = convention.filename(**attributes, startdate='20010101', enddate='20010131')
    # only path
    path     = convention.path(**attributes, startdate='20010101', enddate='20010131')
    # only filename with path
    file     = convention.pattern(root, **attributes, startdate='20010101', enddate='20010131')


Use the conventions module to create your own filenaming conventions:

.. code-block:: python

    from cordex import conventions as conv

    # create your own filename convention string and list
    filename_conv_str  = 'my_convention_{variable}_{model_id}_{domain_id}.nc'
    path_conv_list     = ['model_id','variable']

    # create conventions for filename and path
    filename_conv = conv.FileNameConvention(filename_conv_str)
    path_conv     = conv.FilePathConvention(path_conv_list)


    # now define your attributes to fill the templates.
    root = '/my_root'
    attributes = {'model_id'        : 'GERICS-REMO2015',
                  'variable'        : 'pr',
                  'domain_id'       : 'EUR-11'}

    # create filename and path
    filename = filename_conv.pattern(**attributes)
    path     = path_conv.pattern(root, **attributes)

    # create combined file convention
    file_conv = conv.FileConvention(path_conv, filename_conv)

    # create full filename with path
    file = file_conv.pattern(root, **attributes)

Requirements
------------

* python3.6 or higher
* numpy
* pandas
* (xarray)
* netCDF4
* parse

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
