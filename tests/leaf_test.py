from .context import inlets
import pytest
import os
from shapely.geometry import Point, Polygon
import xarray

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

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
        os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1983-030-0018.che.nc"),
        pytest.param(os.path.join("BOT", "2004-004-0046.che.nc"), marks=pytest.mark.xfail),
    ])
def test_find_temperature_data(source):
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    assert inlets.find_temperature_data(data) is not None

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
        os.path.join("BOT", "1950-050-0040.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        pytest.param(os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail),
    ])
def test_find_salinity_data(source):
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    assert inlets.find_salinity_data(data) is not None

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1978-033-0013.bot.nc"),
        os.path.join("BOT", "1994-022-0001.che.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(os.path.join("BOT", "1976-021-0005.bot.nc"), marks=pytest.mark.xfail),
    ])
def test_find_oxygen_data(source):
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    assert inlets.find_oxygen_data(data) is not None

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
        os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1983-030-0018.che.nc"),
        os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
        os.path.join("BOT", "1950-050-0040.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("BOT", "1978-033-0013.bot.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail),
    ])
def test_find_depth_data(source):
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    assert inlets.find_depth_data(data) is not None

@pytest.mark.parametrize(
    "polygon,point",
    [
        # TODO: Refactor Inlet.contains to take a shapely.Point
        (Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), {"longitude": 0.5, "latitude": 0.5}),
        pytest.param(Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), {"longitude": 2, "latitude": 0.5}, marks=pytest.mark.xfail)
    ])
def test_inlet_contains(polygon, point):
    inlet = inlets.Inlet("Test Inlet", polygon, [1, 2, 3])
    assert inlet.contains(point)

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
        os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1983-030-0018.che.nc"),
        os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
        os.path.join("BOT", "1950-050-0040.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail),
    ])
def test_inlet_add_temperature(source):
    inlet = inlets.Inlet("Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300])
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    inlet.add_temperature_data_from(data)
    assert inlet.has_temperature_data()

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1983-030-0018.che.nc"),
        os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
        os.path.join("BOT", "1950-050-0040.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail),
    ])
def test_inlet_add_salinity(source):
    inlet = inlets.Inlet("Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300])
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    inlet.add_salinity_data_from(data)
    assert inlet.has_salinity_data()

@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1978-033-0013.bot.nc"),
        # TODO: uncomment when umol/kg -> mL/L conversion is complete
        #os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(os.path.join("CTD", "1966-062-0129.ctd.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("CTD", "1966-062-0129.ctd.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("BOT", "1983-030-0018.che.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("BOT", "1950-050-0040.bot.nc"), marks=pytest.mark.xfail),
        pytest.param(os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail),
    ])
def test_inlet_add_oxygen(source):
    inlet = inlets.Inlet("Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300])
    data = xarray.open_dataset(os.path.join(os.path.dirname(__file__), "..", "data", source))
    inlet.add_oxygen_data_from(data)
    assert inlet.has_oxygen_data()
