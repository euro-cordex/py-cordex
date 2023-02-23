"""Unit test package for cordex."""

import importlib

import pytest


def _importorskip(modname):
    try:
        importlib.import_module(modname)
        has = True
    except ImportError:  # pragma: no cover
        has = False
    func = pytest.mark.skipif(not has, reason=f"requires {modname}")
    return has, func


has_matplotlib, requires_matplotlib = _importorskip("matplotlib")
has_cartopy, requires_cartopy = _importorskip("cartopy")
has_regionmask, requires_regionmask = _importorskip("regionmask")
has_xesmf, requires_xesmf = _importorskip("xesmf")
has_geopandas, requires_geopandas = _importorskip("geopandas")
