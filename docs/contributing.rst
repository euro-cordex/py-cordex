.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/euro-cordex/py-cordex/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Cordex Python Package could always use more documentation, whether as part of the
official Cordex Python Package docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/euro-cordex/py-cordex/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `cordex` for local development.

Forking
~~~~~~~

You will need create your own fork of the project. This is really easy using the github
interface, just go to the `py-cordex project page <https://github.com/euro-cordex/py-cordex>`_ and hit the ``Fork`` button.
From your fork, you then clone the repository to your machine::

    git clone https://github.com/your-user-name/py-cordex.git
    cd py-cordex
    git remote add upstream https://github.com/euro-cordex/py-cordex.git


Creating a Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the development environment, we recommend to use the conda package manager.

- Install either `Anaconda <https://www.anaconda.com/download/>`_ or `miniconda
  <https://conda.io/miniconda.html>`_
- Make sure your conda is up to date (``conda update conda``)
- ``cd`` to the *py-cordex* source directory

We don't recommend to use pip installation for development since some
depdencenies (like ``cartopy`` or ``xesmf``) require pre-compiled libraries
in the backend. So the safest way to go is:

1. Install the build dependencies
2. Build and install py-cordex from source

.. code-block:: sh

   # Create and activate the build environment
   conda create -c conda-forge -n py-cordex-tests python=3.9

   conda env update -f ci/requirements/environment.yml

   conda activate py-cordex-tests

   # Build and install py-cordex in editable mode
   pip install -e .

At this point you should be able to import *py-cordex* from your locally
built version:

.. code-block:: sh

   $ python  # start an interpreter
   >>> import cordex
   >>> cordex.__version__

The nice thing about the *editable* mode (that's the ``-e`` flag in the pip install command) is
that you can edit the code directly in the package and use it without having to reinstall
the package. If you work a lot in Jupyter notebooks for development, you should check out
the autoreload magic, e.g., add a cell in the top of your notebook containing:

.. code-block:: sh

   %load_ext autoreload
   %autoreload 2

This will allow you to edit the *py-cordex* source code and use it directly in the notebook
without having to restart the kernel.

See the full conda docs `here <http://conda.pydata.org/docs>`__.
