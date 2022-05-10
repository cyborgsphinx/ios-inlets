from .context import inlets
import numpy
import pytest
import os
from shapely.geometry import Polygon
import xarray


DB_NAME = ":memory:"


def data_path(source):
    return os.path.join(os.path.dirname(__file__), "..", "data", "netCDF_Data", source)


@pytest.mark.parametrize(
    "value,lower,upper",
    [
        (10, 9, 11),
        (11, 9, 11),
        (9, 9, 11),
        pytest.param(8, 9, 11, marks=pytest.mark.xfail),
    ],
)
def test_bounds_check(value, lower, upper):
    assert inlets.is_in_bounds(value, lower, upper)


@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
        os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1983-030-0018.che.nc"),
        pytest.param(
            os.path.join("BOT", "2004-004-0046.che.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_find_temperature_data(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    data = xarray.open_dataset(file)
    assert inlets.find_temperature_data(data) is not None


@pytest.mark.parametrize(
    "source",
    [
        os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
        os.path.join("BOT", "1950-050-0040.bot.nc"),
        os.path.join("CTD", "1966-062-0129.ctd.nc"),
        pytest.param(
            os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_find_salinity_data(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    data = xarray.open_dataset(file)
    assert inlets.find_salinity_data(data) is not None


@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1978-033-0013.bot.nc"),
        os.path.join("BOT", "1994-022-0001.che.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(
            os.path.join("BOT", "1976-021-0005.bot.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_find_oxygen_data(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    data = xarray.open_dataset(file)
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
        os.path.join("BOT", "1994-022-0001.che.nc"),
    ],
)
def test_find_depth_data(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    data = xarray.open_dataset(file)
    assert inlets.find_depth_data(data) is not None


@pytest.mark.parametrize(
    "polygon,point",
    [
        # TODO: Refactor Inlet.contains to take a shapely.Point
        (
            Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]),
            {"longitude": 0.5, "latitude": 0.5},
        ),
        pytest.param(
            Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]),
            {"longitude": 2, "latitude": 0.5},
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_inlet_contains(polygon, point):
    inlet = inlets.Inlet(
        "Test Inlet", "Test Area", polygon, [1, 2, 3], {}, db_name=DB_NAME
    )
    assert inlet.contains(**point)


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
        pytest.param(
            os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_inlet_add_temperature(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    inlet = inlets.Inlet(
        "Test Inlet",
        "Test Area",
        Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]),
        [0, 150, 300],
        {},
        db_name=DB_NAME,
    )
    data = xarray.open_dataset(file)
    inlet.add_data_from_netcdf(data)
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
        os.path.join("BOT", "1994-022-0001.che.nc"),
        pytest.param(
            os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_inlet_add_salinity(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    inlet = inlets.Inlet(
        "Test Inlet",
        "Test Area",
        Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]),
        [0, 150, 300],
        {},
        db_name=DB_NAME,
    )
    data = xarray.open_dataset(file)
    inlet.add_data_from_netcdf(data)
    assert inlet.has_salinity_data()


@pytest.mark.parametrize(
    "source",
    [
        os.path.join("BOT", "1930-031-0001.bot.nc"),
        os.path.join("CTD", "2021-020-0001.ctd.nc"),
        os.path.join("BOT", "1978-033-0013.bot.nc"),
        os.path.join("BOT", "1994-022-0001.che.nc"),
        os.path.join("BOT", "1994-003-0001.che.nc"),
        pytest.param(
            os.path.join("CTD", "1966-062-0129.ctd.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("CTD", "1966-062-0129.ctd.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
            marks=pytest.mark.xfail,
        ),
        pytest.param(
            os.path.join("CUR", "CM1_19890407_19890504_0020m.cur.nc"),
            marks=pytest.mark.xfail,
        ),
        pytest.param(
            os.path.join("BOT", "1983-030-0018.che.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("CUR", "C_19780823_19781204_0175m.cur.nc"),
            marks=pytest.mark.xfail,
        ),
        pytest.param(
            os.path.join("BOT", "1950-050-0040.bot.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_inlet_add_oxygen(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    inlet = inlets.Inlet(
        "Test Inlet",
        "Test Area",
        Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]),
        [0, 150, 300],
        {},
        db_name=DB_NAME,
    )
    data = xarray.open_dataset(file)
    inlet.add_data_from_netcdf(data)
    assert inlet.has_oxygen_data()


@pytest.mark.parametrize(
    "data,placeholder",
    [
        ([0.0, 1.0, -99.0], -99),
        pytest.param([0.0, 1.0, 2.0], -99, marks=pytest.mark.xfail),
    ],
)
def test_reinsert_nan(data, placeholder):
    assert any(numpy.isnan(inlets.reinsert_nan(data, placeholder)))
