import argparse
import datetime
import inlets
import itertools
import matplotlib.pyplot as plt
import os
from typing import List

END = datetime.datetime.now()
INLET_LINES = ["m-s", "y-d", "k-o", "c-^", "b-d", "g-s", "r-s"]


###################
# Utility functions
###################

def normalize(string: str):
    return string.strip().lower().replace(" ", "-")


def figure_path(filename: str):
    return os.path.join("figures", filename)


########################
# Single inlet functions
########################


def chart_data(inlet: inlets.Inlet, limits: List[float], data_fn):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    shallow_time, shallow_data = data_fn(inlet, inlets.SHALLOW)
    middle_time, middle_data = data_fn(inlet, inlets.MIDDLE)
    deep_time, deep_data = data_fn(inlet, inlets.DEEP)
    if len(limits) > 1:
        shallow_time, shallow_data = zip(
            *[
                [t, d]
                for t, d in zip(shallow_time, shallow_data)
                if limits[0] < d < limits[1]
            ]
        )
        middle_time, middle_data = zip(
            *[
                [t, d]
                for t, d in zip(middle_time, middle_data)
                if limits[0] < d < limits[1]
            ]
        )
        deep_time, deep_data = zip(
            *[[t, d] for t, d in zip(deep_time, deep_data) if limits[0] < d < limits[1]]
        )
    plt.plot(
        shallow_time,
        shallow_data,
        "xg",
        label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m",
    )
    plt.plot(
        middle_time,
        middle_data,
        "+m",
        label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m",
    )
    deep_label = (
        f">{inlet.deep_bounds[0]}m"
        if inlet.deep_bounds[1] is None
        else f"{inlet.deep_bounds[0]}m-{inlet.deep_bounds[1]}m"
    )
    plt.plot(deep_time, deep_data, "xb", label=deep_label)
    plt.legend()


def chart_temperatures(inlet: inlets.Inlet, limits: List[float], use_averages: bool):
    chart_data(
        inlet,
        limits,
        lambda inlet, bucket: inlet.get_temperature_data(
            bucket, before=END, do_average=use_averages
        ),
    )
    plt.ylabel("Temperature (C)")
    plt.title(f"{inlet.name} Deep Water Temperatures")


def chart_salinities(inlet: inlets.Inlet, limits: List[float], use_averages: bool):
    chart_data(
        inlet,
        limits,
        lambda inlet, bucket: inlet.get_salinity_data(
            bucket, before=END, do_average=use_averages
        ),
    )
    plt.ylabel("Salinity (PSU)")
    plt.title(f"{inlet.name} Deep Water Salinity")


def chart_oxygen_data(inlet: inlets.Inlet, limits: List[float], use_averages: bool):
    chart_data(
        inlet,
        limits,
        lambda inlet, bucket: inlet.get_oxygen_data(
            bucket, before=END, do_average=use_averages
        ),
    )
    plt.ylabel("DO (ml/l)")
    plt.title(f"{inlet.name} Deep Water Dissolved Oxygen")


def chart_stations(inlet: inlets.Inlet, _limits: List[float], _use_averages: bool):
    # `limits` and `use_averages` are only present to conform to expected function type
    del _limits, _use_averages
    plt.clf()
    data = []
    for year, stations in inlet.get_station_data(before=END).items():
        data.extend([year for _ in range(len(stations))])
    n_bins = max(data) - min(data) + 1
    plt.hist(data, bins=n_bins, align="left", label=f"Number of files {len(data)}")
    plt.ylabel("Number of Stations")
    plt.legend()
    plt.title(f"{inlet.name} Sampling History")


def do_chart(
    inlet: inlets.Inlet,
    kind: str,
    use_limits: bool,
    chart_fn,
    use_averages: bool,
):
    print(f"Producing {kind} plot for {inlet.name}")
    chart_fn(
        inlet,
        inlet.limits[kind] if use_limits and kind in inlet.limits else [],
        use_averages,
    )
    limits = "" if use_limits else "-full"
    average = "-average" if use_averages else ""
    plt.savefig(figure_path(f"{normalize(inlet.name)}-{kind}{limits}{average}.png"))

#####################
# All inlet functions
#####################

def chart_all_data(times, data, label=""):
    plt.scatter(times, data, label=label)


def bounds_label(inlet, bucket):
    bounds_name = bucket.lower() + "_bounds"
    bounds = getattr(inlet, bounds_name)
    label = inlet.name
    if bounds[1] is None:
        label += f" >{bounds[0]}"
    else:
        label += f" {bounds[0]}-{bounds[1]}"
    return label


def chart_all_temperature(inlet, bucket):
    label = bounds_label(inlet, bucket)
    chart_all_data(
        *inlet.get_temperature_data(bucket, before=END, do_average=True), label=label
    )
    plt.ylabel("Temperature (C)")


def chart_all_salinity(inlet, bucket):
    label = bounds_label(inlet, bucket)
    chart_all_data(
        *inlet.get_salinity_data(bucket, before=END, do_average=True), label=label
    )
    plt.ylabel("Salinity (PSU)")


def chart_all_oxygen(inlet, bucket):
    label = bounds_label(inlet, bucket)
    chart_all_data(
        *inlet.get_oxygen_data(bucket, before=END, do_average=True), label=label
    )
    plt.ylabel("DO (ml/l)")


def do_chart_all(inlet_list, kind, bucket, chart_all_fn):
    print(f"Producing {kind} plot for {bucket}")
    plt.clf()
    for inlet in inlet_list:
        chart_all_fn(inlet, bucket)
    names = "-".join(normalize(inlet.name) for inlet in inlet_list)
    bounds_name = bucket.lower() + "_bounds"
    lowest = min(getattr(inlet, bounds_name)[0] for inlet in inlet_list)
    highest = max(
        getattr(inlet, bounds_name)[1]
        if getattr(inlet, bounds_name)[1] is not None
        else 0
        for inlet in inlet_list
    )
    if highest < lowest:
        plt.title(f"{kind.capitalize()} comparison below {lowest}m")
        bounds = f"{lowest}-bottom"
    else:
        plt.title(f"{kind.capitalize()} comparison from {lowest}m to {highest}m")
        bounds = f"{lowest}-{highest}"
    plt.legend()
    plt.savefig(figure_path(f"{bounds}-{kind}-{names}.png"))

############################
# Annual averaging functions
############################

def do_annual_work(inlet_list, data_fn, averaging_fn):
    plt.clf()
    for inlet, line_style in zip(inlet_list, INLET_LINES):
        totals = {}
        times, data = data_fn(inlet)
        for time, datum in zip(times, data):
            year = time.year
            if year not in totals:
                totals[year] = (0, 0)
            total, num = totals[year]
            totals[year] = (total + datum, num + 1)

        averages = averaging_fn(totals, data)
        years, values = zip(*sorted(averages.items(), key=lambda item: item[0]))
        plt.plot(years, values, line_style, label=inlet.name)

    plt.legend()


def do_annual_work_parts(inlet_list, data_fn, averaging_fn, y_label, title):
    offset = 0
    loops = 0
    stride = len(INLET_LINES)
    while offset < len(inlet_list):
        loops += 1
        do_annual_work(inlet_list[offset:offset+stride], data_fn, averaging_fn)
        offset += stride
        plt.ylabel(y_label)
        plt.title(title)
        plt.savefig(figure_path(f"{normalize(title)}-part{loops}.png"))


def annual_averaging(totals, _data):
    del _data
    return {y: t / n for y, (t, n) in totals.items()}


def anomalies_averaging(totals, data):
    avg = sum(data) / len(data)
    avgs = annual_averaging(totals, data)
    return {y: a - avg for y, a in avgs.items()}


def chart_anomalies(inlet_list: List[inlets.Inlet], data_fn, y_label, title):
    do_annual_work_parts(inlet_list, data_fn, anomalies_averaging, y_label, title)


def chart_temperature_anomalies(inlet_list: List[inlets.Inlet]):
    print("Producing temperature anomaly plots")
    chart_anomalies(inlet_list, lambda inlet: inlet.get_temperature_data(inlets.ALL), "Temperature (C)", "Deep Water Temperature Anomalies")


def chart_salinity_anomalies(inlet_list: List[inlets.Inlet]):
    print("Producing salinity anomaly plots")
    chart_anomalies(inlet_list, lambda inlet: inlet.get_salinity_data(inlets.ALL), "Salinity (PSU)", "Deep Water Salinity Anomalies")


def chart_oxygen_anomalies(inlet_list: List[inlets.Inlet]):
    print("Producing oxygen anomaly plot")
    chart_anomalies(inlet_list, lambda inlet: inlet.get_oxygen_data(inlets.ALL), "Oxygen (mL/L)", "Deep Water Dissolved Oxygen Anomalies")


def do_chart_annual_averages(inlet_list: List[inlets.Inlet], data_fn, y_label, title):
    do_annual_work_parts(inlet_list, data_fn, annual_averaging, y_label, title)


def chart_annual_temperature_averages(inlet_list: List[inlets.Inlet]):
    print("Producing annual temperature plots")
    do_chart_annual_averages(inlet_list, lambda inlet: inlet.get_temperature_data(inlets.ALL), "Temperature (C)", "Deep Water Temperature Annual Averages")


def chart_annual_salinity_averages(inlet_list: List[inlets.Inlet]):
    print("Producing annual salinity plots")
    do_chart_annual_averages(inlet_list, lambda inlet: inlet.get_salinity_data(inlets.ALL), "Salinity (PSU)", "Deep Water Salinity Annual Averages")


def chart_annual_oxygen_averages(inlet_list: List[inlets.Inlet]):
    print("Producing annual oxygen plots")
    do_chart_annual_averages(inlet_list, lambda inlet: inlet.get_oxygen_data(inlets.ALL), "Oxygen (mL/L)", "Deep Water Dissolved Oxygen Annual Averages")


###################
# Monthly functions
###################

def chart_monthly_sample(inlet: inlets.Inlet):
    months = [
        "padding",
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ]
    files = {
        "January": {},
        "February": {},
        "March": {},
        "April": {},
        "May": {},
        "June": {},
        "July": {},
        "August": {},
        "September": {},
        "October": {},
        "November": {},
        "December": {},
    }
    min_year = END.year
    max_year = 0
    for datum in itertools.chain(
        inlet.data.get_temperature_data(),
        inlet.data.get_salinity_data(),
        inlet.data.get_oxygen_data(),
    ):
        year = datum.time.year
        min_year = min(year, min_year)
        max_year = max(year, max_year)
        month = months[datum.time.month]
        if year not in files[month]:
            files[month][year] = set()
        files[month][year].add(datum.filename)
    year_range = max_year - min_year + 1
    values = []
    for _ in range(year_range):
        values.append([0] * 12)
    for month, d in files.items():
        for year, filenames in d.items():
            year_idx = year - min_year
            month_idx = months.index(month) - 1
            values[year_idx][month_idx] += len(filenames)
    biggest = 0
    for row in values:
        biggest = max(biggest, *row)

    plt.clf()
    plt.figure(figsize=(40, 10), constrained_layout=True)
    plt.imshow(list(map(list, zip(*values))), vmin=0, vmax=biggest, cmap="Blues")
    plt.yticks(ticks=range(12), labels=months[1:])
    plt.xticks(ticks=range(year_range), labels=range(min_year, max_year + 1))
    for i, _ in enumerate(values):
        for j, _ in enumerate(values[i]):
            plt.text(i, j, values[i][j], ha="center", va="center", color="k")

    plt.title(f"{inlet.name} Sampling frequency by month")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-monthly-sampling.png"))

def main():
    parser = argparse.ArgumentParser()
    # inlet retrieval args
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-n", "--skip-netcdf", action="store_true")
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    # plot args
    parser.add_argument("-l", "--no-limits", action="store_true")
    parser.add_argument("-i", "--inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-k", "--limit-name", type=str, nargs="+", default=[])
    parser.add_argument("-I", "--remove-inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-b", "--plot-buckets", action="store_true")
    parser.add_argument("-a", "--use-averages", action="store_true")
    parser.add_argument("-A", "--plot-annual", action="store_true")
    parser.add_argument("-m", "--plot-monthly", action="store_true")
    args = parser.parse_args()
    inlet_list = inlets.get_inlets(
        args.data,
        args.from_saved,
        args.skip_netcdf,
        args.inlet_name,
        args.remove_inlet_name,
        args.limit_name,
    )
    plt.figure(figsize=(8, 6))
    if args.plot_annual:
        chart_annual_temperature_averages(inlet_list)
        chart_annual_salinity_averages(inlet_list)
        chart_annual_oxygen_averages(inlet_list)
        chart_temperature_anomalies(inlet_list)
        chart_salinity_anomalies(inlet_list)
        chart_oxygen_anomalies(inlet_list)
    elif args.plot_monthly:
        for inlet in inlet_list:
            chart_monthly_sample(inlet)
    elif args.plot_buckets:
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.SHALLOW,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.MIDDLE,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.DEEP,
            chart_all_temperature,
        )
        do_chart_all(inlet_list, "salinity", inlets.SHALLOW, chart_all_salinity)
        do_chart_all(inlet_list, "salinity", inlets.MIDDLE, chart_all_salinity)
        do_chart_all(inlet_list, "salinity", inlets.DEEP, chart_all_salinity)
        do_chart_all(inlet_list, "oxygen", inlets.SHALLOW, chart_all_oxygen)
        do_chart_all(inlet_list, "oxygen", inlets.MIDDLE, chart_all_oxygen)
        do_chart_all(inlet_list, "oxygen", inlets.DEEP, chart_all_oxygen)
    else:
        for inlet in inlet_list:
            do_chart(
                inlet,
                "temperature",
                not args.no_limits,
                chart_temperatures,
                args.use_averages,
            )
            do_chart(
                inlet,
                "salinity",
                not args.no_limits,
                chart_salinities,
                args.use_averages,
            )
            do_chart(
                inlet,
                "oxygen",
                not args.no_limits,
                chart_oxygen_data,
                args.use_averages,
            )
            do_chart(
                inlet,
                "samples",
                not args.no_limits,
                chart_stations,
                args.use_averages,
            )
    plt.close()


if __name__ == "__main__":
    main()
