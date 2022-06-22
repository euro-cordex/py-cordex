#!/usr/bin/env python
from setuptools import setup

# get version
# with open("cordex/version.py") as f:
#    line = f.readline().strip().replace(" ", "").replace('"', "")
#    version = line.split("=")[1]
#    __version__ = version
#
# setup(version=__version__)
setup(use_scm_version={"fallback_version": "999"})
