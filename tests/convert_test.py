#!/usr/bin/env python3

from .context import inlets
import pytest

import convert


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
    actuals = convert.convert_umol_kg_to_mL_L(
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
    actuals = convert.convert_umol_kg_to_mL_L(oxygen_umol_kg, longitude, latitude)
    assert len(actuals[0]) == len(expecteds)
    for actual, expected in zip(actuals[0], expecteds):
        assert abs(actual - expected) < 0.21


@pytest.mark.parametrize(
    "oxygen_percent,temperature_C,salinity_SP,expecteds",
    [
        (
            [138.01, 140.83, 92.70, 71.14, 53.76, 52.05, 51.67, 50.75],
            [10.1993, 9.9926, 8.6750, 8.6019, 8.0816, 8.0756, 8.0672, 8.0458],
            [30.2153, 30.7551, 31.9281, 32.0355, 32.8430, 32.8635, 32.8692, 32.8796],
            # above values given to https://ocena.ices.dk/tools/calculator.aspx to produce expected values
            [8.949, 9.143, 6.152, 4.726, 3.595, 3.481, 3.456, 3.396],
        )
    ],
)
def test_convert_oxygen_percent_to_ml_per_l(
    oxygen_percent, temperature_C, salinity_SP, expecteds
):
    actuals = convert.convert_percent_to_mL_L(
        oxygen_percent, temperature_C, salinity_SP
    )
    assert actuals is not None
    assert len(actuals) == len(expecteds)
    for actual, expected in zip(actuals, expecteds):
        assert abs(actual - expected) < 0.001
