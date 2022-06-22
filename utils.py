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
    return math.sqrt(sum(map(lambda x: (x - u) ** 2, data)) / len(data))


def index_by_month(dates):
    dates = list(dates)
    start_year = min(date.year for date in dates)
    return [date.month + (12 * (date.year - start_year)) for date in dates]


def remove_seasonal_trend(x, y, remove_trend=False, by_difference=True, remove_sd=True):
    to_process = list(y)
    if remove_trend:
        if by_difference:
            to_process = list(map(lambda a, b: a - b, to_process[1:], to_process))
        else:
            domain = index_by_month(x)
            fit = Polynomial.fit(domain, to_process, 1)
            coeficients = fit.convert().coef
            to_process = list(
                map(
                    lambda a, b: a - (coeficients[0] + coeficients[1] * b),
                    to_process,
                    domain,
                )
            )
    avg = mean(to_process)
    out = [x - avg for x in to_process]
    if remove_sd:
        std_dev = sd(to_process)
        out = [x / std_dev for x in out]
    return out
