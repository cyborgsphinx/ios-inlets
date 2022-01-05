from .context import inlets
import numpy
import pytest
import os
from shapely.geometry import Polygon
import xarray


def data_path(source):
    return os.path.join(os.path.dirname(__file__), "..", "data", "netCDF_Data", source)


@pytest.mark.parametrize(
    "value,lower,upper",
    [
        (10, 9, 11),
        (11, 9, 11),
        pytest.param(9, 9, 11, marks=pytest.mark.xfail),
    ],
)
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
    ],
)
def test_find_any(source, attrs, expected):
    assert inlets.find_any(source, *attrs) == expected


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
        pytest.param(
            os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail
        ),
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
        "Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300]
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
        pytest.param(
            os.path.join("BOT", "1978-033-0013.bot.nc"), marks=pytest.mark.xfail
        ),
        pytest.param(
            os.path.join("ADCP", "nep1_20060512_20060525_0095m.adcp.L1.nc"),
            marks=pytest.mark.xfail,
        ),
        pytest.param(
            os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_inlet_add_salinity(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    inlet = inlets.Inlet(
        "Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300]
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
        # TODO: uncomment when umol/kg -> mL/L conversion is complete
        # os.path.join("BOT", "1994-003-0001.che.nc"),
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
        pytest.param(
            os.path.join("BOT", "1994-022-0001.che.nc"), marks=pytest.mark.xfail
        ),
    ],
)
def test_inlet_add_oxygen(source):
    file = data_path(source)
    if not os.path.isfile(file):
        pytest.skip(f"{file} is not accessible")
    inlet = inlets.Inlet(
        "Test Inlet", Polygon([[0, 0], [0, 1], [1, 1], [1, 0]]), [0, 150, 300]
    )
    data = xarray.open_dataset(file)
    inlet.add_data_from_netcdf(data)
    assert inlet.has_oxygen_data()


@pytest.mark.parametrize(
    "oxygen_umol_kg,temperature_C,salinity_SP,pressure_dbar,longitude,latitude,expecteds",
    [
        (
            [
                201.3,
                192.5,
                172.2,
                161.3,
                145.7,
                115.6,
                139.9,
                138.9,
                120.9,
                101.2,
                75.5,
                65.1,
                61.2,
                60.0,
            ],
            [
                11.0228,
                10.6766,
                9.827,
                9.4797,
                9.2095,
                8.1733,
                7.5017,
                7.3041,
                7.1969,
                6.9334,
                6.9748,
                6.8451,
                6.7665,
                6.7248,
            ],
            [
                31.7055,
                31.7774,
                31.9284,
                32.1258,
                32.1665,
                32.8476,
                33.12,
                33.4193,
                33.5858,
                33.7288,
                33.8478,
                33.8913,
                33.9206,
                33.9263,
            ],
            [
                1.1,
                5.1,
                9.9,
                19.9,
                30.4,
                40.8,
                50.8,
                75.3,
                99.9,
                125.0,
                150.2,
                175.3,
                200.7,
                249.0,
            ],
            -124.73517,
            48.50083,
            [
                4.62,
                4.42,
                3.95,
                3.70,
                3.34,
                2.66,
                3.21,
                3.19,
                2.78,
                2.33,
                1.74,
                1.50,
                1.41,
                1.38,
            ],
        ),
    ],
)
def test_convert_oxygen_data_full_data(
    oxygen_umol_kg,
    temperature_C,
    salinity_SP,
    pressure_dbar,
    longitude,
    latitude,
    expecteds,
):
    actuals = inlets.convert_umol_kg_to_mL_L(
        oxygen_umol_kg, longitude, latitude, temperature_C, salinity_SP, pressure_dbar
    )
    assert len(actuals[0]) == len(expecteds)
    for actual, expected, oxygen, temperature, salinity, pressure in zip(
        actuals[0], expecteds, oxygen_umol_kg, temperature_C, salinity_SP, pressure_dbar
    ):
        assert (
            abs(actual - expected) < 0.01
        ), f"Failure in conversion at oxygen(umol/kg)={oxygen}, temperature(C)={temperature}, salinity(PSU)={salinity}, pressure(dbar)={pressure}"


@pytest.mark.parametrize(
    "oxygen_umol_kg,longitude,latitude,expecteds",
    [
        (
            [291.6, 292.8, 240.9, 194.6, 184.2, 172.8, 115.0, 110.1],
            -124.2405,
            49.7585,
            [6.73, 6.711, 5.508, 4.455, 4.218, 3.957, 2.635, 2.523],
        )
    ],
)
def test_convert_oxygen_data_missing_data(
    oxygen_umol_kg, longitude, latitude, expecteds
):
    actuals = inlets.convert_umol_kg_to_mL_L(oxygen_umol_kg, longitude, latitude)
    assert len(actuals[0]) == len(expecteds)
    for actual, expected in zip(actuals[0], expecteds):
        assert abs(actual - expected) < 0.21


@pytest.mark.parametrize(
    "data,placeholder",
    [
        ([0.0, 1.0, -99.0], -99),
        pytest.param([0.0, 1.0, 2.0], -99, marks=pytest.mark.xfail),
    ],
)
def test_reinsert_nan(data, placeholder):
    assert any(numpy.isnan(inlets.reinsert_nan(data, placeholder)))


@pytest.mark.parametrize(
    "source,prefix,index",
    [
        (["depth", "temperature"], "depth", 0),
        (["depth", "temperature"], "temperature", 1),
        (["depth", "temperature"], "oxygen", -1),
        (["depth:suffix", "temperature"], "depth", 0),
    ],
)
def test_find_first(source, prefix, index):
    assert inlets.find_first(source, prefix) == index
