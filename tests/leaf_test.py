from .context import inlets
import pytest
from shapely.geometry import Point, Polygon

@pytest.mark.parametrize(
    "name,lat,long",
    [
        ("Saanich Inlet",  -123.492, 48.6583),
        ("Jervis Inlet",   -124.043, 49.8017),
        ("Bute Inlet",     -124.969, 50.5318),
        ("Knight Inlet",   -125.781, 50.6973),
        ("Indian Arm",     -123.087, 49.2987),
        ("Howe Sound",     -123.294, 49.4779),
        ("Muchalat Inlet", -126.176, 49.6639),
    ])
def test_find_polygon(name, lat, long):
    polygon = inlets.find_polygon_for(name)
    assert Polygon(polygon).contains(Point(lat, long))

@pytest.mark.parametrize(
    "value,lower,upper",
    [
        (10, 9, 11),
        (11, 9, 11),
        pytest.param(9, 9, 11, marks=pytest.mark.xfail),
    ])
def test_bounds_check(value, lower, upper):
    assert inlets.is_in_bounds(value, lower, upper)

class Data(object):
    def __init__(self, val):
        self.val = val

@pytest.mark.parametrize(
    "source,attrs,expected",
    [
        (Data(9), ["val"], 9),
        (Data(9), ["foo"], None),
        (Data(9), ["foo", "val"], 9),
        (Data(9), ["val", "bar"], 9),
        (Data(9), ["foo", "bar"], None),
    ])
def test_find_any(source, attrs, expected):
    assert inlets.find_any(source, attrs) == expected
