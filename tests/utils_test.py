#!/usr/bin/env python3

from .context import inlets
import pytest
import datetime

import utils


@pytest.mark.parametrize(
    "arg,expected",
    [
        ("Saanich Inlet", "saanich-inlet"),
        ("saanich-inlet", "saanich-inlet"),
    ],
)
def test_normalize(arg, expected):
    assert utils.normalize(arg) == expected


@pytest.mark.parametrize(
    "lower,upper,expected",
    [
        (0, 30, "0m-30m"),
        (100, None, ">100m"),
    ],
)
def test_label_from_bounds(lower, upper, expected):
    assert utils.label_from_bounds(lower, upper) == expected


@pytest.mark.parametrize(
    "totals,key,datum,expected",
    [
        ({}, 1, 2, (2, 1)),
        ({2: "other"}, 1, 2, (2, 1)),
        ({1: (3, 1)}, 1, 2, (5, 2)),
    ],
)
def test_update_totals(totals, key, datum, expected):
    utils.update_totals(totals, key, datum)
    assert totals[key] == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 5.5),
    ],
)
def test_mean(data, expected):
    assert utils.mean(data) == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        ([1, 2, 1, 2, 1, 2, 1, 2, 1, 2], 0.5),
    ],
)
def test_sd(data, expected):
    assert utils.sd(data) == expected


@pytest.mark.parametrize(
    "dates,expected",
    [
        (
            [
                datetime.date(1950, 1, 1),
                datetime.date(1950, 2, 1),
                datetime.date(1951, 1, 1),
            ],
            [1, 2, 13],
        ),
    ],
)
def test_index_by_month(dates, expected):
    assert utils.index_by_month(dates) == expected
