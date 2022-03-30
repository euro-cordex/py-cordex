import pytest
from cordex import cmor
import cftime as cfdt
import datetime as dt


def test_cftime():
    assert cmor.to_cftime(dt.datetime(2000, 1, 1, 1)) == cfdt.datetime(2000, 1, 1, 1)


@pytest.mark.parametrize("dt", [dt, cfdt])
def test_season(dt):
    assert cmor.season(dt.datetime(2000, 1, 1)) == "DJF"
    assert cmor.season(dt.datetime(2001, 3, 1)) == "MAM"
    assert cmor.season(dt.datetime(2001, 6, 1)) == "JJA"
    assert cmor.season(dt.datetime(2001, 9, 1)) == "SON"
    assert cmor.season(dt.datetime(2001, 12, 1)) == "DJF"
    bounds = (dt.datetime(1999, 12, 1, 0, 0), dt.datetime(2000, 3, 1, 0, 0))
    assert cmor.season_bounds(dt.datetime(2000, 1, 1)) == bounds
    assert cmor.mid_of_season(dt.datetime(2000, 1, 1)) == dt.datetime(
        2000, 1, 15, 12, 0
    )


def test_month():
    assert cmor.month_bounds(dt.datetime(2000, 1, 10)) == (
        dt.datetime(2000, 1, 1),
        dt.datetime(2000, 2, 1),
    )
    assert cmor.month_bounds(dt.datetime(2000, 12, 10)) == (
        dt.datetime(2000, 12, 1),
        dt.datetime(2001, 1, 1),
    )
    assert cmor.mid_of_month(dt.datetime(2000, 1, 1)) == dt.datetime(2000, 1, 16, 12)
    # leap years
    assert cmor.mid_of_month(dt.datetime(2000, 2, 1)) == dt.datetime(2000, 2, 15, 12)
    assert cmor.mid_of_month(dt.datetime(2001, 2, 1)) == dt.datetime(2001, 2, 15)


def test_cfmonth():
    assert cmor.month_bounds(dt.datetime(2000, 1, 10)) == (
        dt.datetime(2000, 1, 1),
        dt.datetime(2000, 2, 1),
    )
    # leap years
    assert cmor.mid_of_month(cfdt.datetime(2000, 2, 1)) == cfdt.datetime(
        2000, 2, 15, 12
    )
    assert cmor.mid_of_month(cfdt.datetime(2001, 2, 1)) == cfdt.datetime(2001, 2, 15)
    assert cmor.mid_of_month(
        cfdt.datetime(2000, 2, 1, calendar="360_day")
    ) == cfdt.datetime(2000, 2, 16, calendar="360_day")
    assert cmor.mid_of_month(
        cfdt.datetime(2001, 2, 1, calendar="360_day")
    ) == cfdt.datetime(2001, 2, 16, calendar="360_day")
