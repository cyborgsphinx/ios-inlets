import datetime
from enum import Enum
import math
from numpy.polynomial.polynomial import Polynomial


def normalize(string: str):
    return string.strip().lower().replace(" ", "-")


def label_from_bounds(lower, upper):
    if upper is None:
        return f">{lower}m"
    else:
        return f"{lower}m-{upper}m"


def update_totals(totals, key, datum):
    if key not in totals:
        totals[key] = (0, 0)
    total, num = totals[key]
    totals[key] = (total + datum, num + 1)


def mean(data):
    return sum(data) / len(data)


def sd(data):
    u = mean(data)
    return math.sqrt(sum((x - u) ** 2 for x in data) / len(data))


def index_by_month(dates):
    dates = list(dates)
    start_year = min(date.year for date in dates)
    return [date.month + (12 * (date.year - start_year)) for date in dates]


def date_from_float(num):
    # ignoring leap years for simplicity
    year = math.trunc(num)
    # want 1-365 instead of 0-364
    day = math.trunc((num - year) * 365) + 1
    months = [
        0,  # padding
        31,  # January
        28,  # February
        31,  # March
        30,  # April
        31,  # May
        30,  # June
        31,  # July
        31,  # August
        30,  # September
        31,  # October
        30,  # November
        31,  # December
    ]
    for m, d in enumerate(months):
        if day < d:
            month = m
            break
        else:
            day -= d
    return datetime.date(year, month, day)


class Trend(Enum):
    NONE = 0
    DIFF = 1
    LINEAR = 2


def remove_seasonal_trend(x, y, trend=Trend.NONE, remove_sd=True):
    to_process = list(y)
    if trend == Trend.DIFF:
        to_process = [a - b for a, b in zip(to_process[1:], to_process)]
    elif trend == Trend.LINEAR:
        domain = index_by_month(x)
        fit = Polynomial.fit(domain, to_process, 1)
        coeficients = fit.convert().coef
        to_process = [
            a - (coeficients[0] + coeficients[1] * b)
            for a, b in zip(to_process, domain)
        ]
    avg = mean(to_process)
    out = [x - avg for x in to_process]
    if remove_sd:
        std_dev = sd(to_process)
        out = [x / std_dev for x in out]
    return out
