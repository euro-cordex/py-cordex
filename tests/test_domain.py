# -*- coding: utf-8 -*-
# flake8: noqa
import pytest
from cordex import cordex_domain as dm


def test_write():
    domain = dm('EUR-11')
    domain.to_netcdf('EUR-11.nc')
