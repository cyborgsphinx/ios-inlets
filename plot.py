import argparse
import datetime
import inlets
import matplotlib.pyplot as plt
import numpy
import os
from typing import List

END = datetime.datetime.now()
INLET_LINES = ["m-s", "y-d", "k-o", "c-^", "b-d", "g-s", "r-s"]


def chart_data(inlet: inlets.Inlet, limits: List[float], data_fn):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    if len(limits) > 1:
        axes = plt.axes()
        axes.set_ylim(limits[0], limits[1])
    shallow_time, shallow_data = data_fn(inlet, inlets.SHALLOW)
    plt.plot(
        shallow_time,
        shallow_data,
        "xg",
        label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m",
    )
    middle_time, middle_data = data_fn(inlet, inlets.MIDDLE)
    plt.plot(
        middle_time,
        middle_data,
        "+m",
        label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m",
    )
    deep_time, deep_data = data_fn(inlet, inlets.DEEP)
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
    plt.clf()
    data = []
    for year, stations in inlet.get_station_data(before=END).items():
        data.extend([year for _ in range(len(stations))])
    n_bins = max(data) - min(data) + 1
    plt.hist(data, bins=n_bins, align="left", label=f"Number of files {len(data)}")
    plt.ylabel("Number of Stations")
    plt.legend()
    plt.title(f"{inlet.name} Sampling History")


def normalize(string: str):
    return string.strip().lower().replace(" ", "-")


def do_chart(
    inlet: inlets.Inlet,
    kind: str,
    show_figure: bool,
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
    if show_figure:
        plt.show()
    else:
        limits = "" if use_limits else "-full"
        plt.savefig(
            os.path.join("figures", f"{normalize(inlet.name)}-{kind}{limits}.png")
        )


def chart_anomalies(inlet_list: List[inlets.Inlet], data_fn):
    plt.clf()
    for inlet, line_style in zip(inlet_list, INLET_LINES):
        totals = {}
        data = [
            datum
            for datum in data_fn(inlet)
            if datum.bucket != inlets.IGNORE and not numpy.isnan(datum.datum)
        ]
        for datum in data:
            year = datum.time.year
            if year not in totals:
                totals[year] = (0, 0)
            total, num = totals[year]
            totals[year] = (total + datum.datum, num + 1)

        avg = sum(x.datum for x in data) / len(data)
        avgs = {y: t / n for y, (t, n) in totals.items()}

        avg_diffs = {y: a - avg for y, a in avgs.items()}

        years, anomalies = zip(*sorted(avg_diffs.items(), key=lambda item: item[0]))
        plt.plot(years, anomalies, line_style, label=inlet.name)

    plt.legend()


def chart_temperature_anomalies(inlet_list: List[inlets.Inlet], show_figure: bool):
    print("Producing temperature anomaly plot")
    chart_anomalies(inlet_list, lambda inlet: inlet.temperature_data)

    plt.ylabel("Temperature (C)")
    plt.title("Deep Water Temperature Anomalies")
    if show_figure:
        plt.show()
    else:
        plt.savefig(os.path.join("figures", f"deep-water-temperature-anomalies.png"))


def chart_salinity_anomalies(inlet_list: List[inlets.Inlet], show_figure: bool):
    print("Producing salinity anomaly plot")
    chart_anomalies(inlet_list, lambda inlet: inlet.salinity_data)

    plt.ylabel("Salinity (PSU)")
    plt.title("Deep Water Salinity Anomalies")
    if show_figure:
        plt.show()
    else:
        plt.savefig(os.path.join("figures", f"deep-water-salinity-anomalies.png"))


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


def do_chart_all(inlet_list, kind, bucket, show_figure, chart_all_fn):
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
    if show_figure:
        plt.show()
    else:
        plt.savefig(os.path.join("figures", f"{bounds}-{kind}-{names}.png"))


def main():
    parser = argparse.ArgumentParser()
    # inlet retrieval args
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-n", "--skip-netcdf", action="store_true")
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    # plot args
    parser.add_argument("-s", "--show-figure", action="store_true")
    parser.add_argument("-l", "--no-limits", action="store_true")
    parser.add_argument("-i", "--inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-k", "--limit-name", type=str, nargs="+", default=[])
    parser.add_argument("-I", "--remove-inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-b", "--plot-buckets", action="store_true")
    parser.add_argument("-a", "--use-averages", action="store_true")
    parser.add_argument("-p", "--plot-anomalies", action="store_true")
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
    if args.plot_buckets:
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.SHALLOW,
            args.show_figure,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.MIDDLE,
            args.show_figure,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.DEEP,
            args.show_figure,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list, "salinity", inlets.SHALLOW, args.show_figure, chart_all_salinity
        )
        do_chart_all(
            inlet_list, "salinity", inlets.MIDDLE, args.show_figure, chart_all_salinity
        )
        do_chart_all(
            inlet_list, "salinity", inlets.DEEP, args.show_figure, chart_all_salinity
        )
        do_chart_all(
            inlet_list, "oxygen", inlets.SHALLOW, args.show_figure, chart_all_oxygen
        )
        do_chart_all(
            inlet_list, "oxygen", inlets.MIDDLE, args.show_figure, chart_all_oxygen
        )
        do_chart_all(
            inlet_list, "oxygen", inlets.DEEP, args.show_figure, chart_all_oxygen
        )
    else:
        for inlet in inlet_list:
            do_chart(
                inlet,
                "temperature",
                args.show_figure,
                not args.no_limits,
                chart_temperatures,
                args.use_averages,
            )
            do_chart(
                inlet,
                "salinity",
                args.show_figure,
                not args.no_limits,
                chart_salinities,
                args.use_averages,
            )
            do_chart(
                inlet,
                "oxygen",
                args.show_figure,
                not args.no_limits,
                chart_oxygen_data,
                args.use_averages,
            )
            do_chart(
                inlet,
                "stations",
                args.show_figure,
                not args.no_limits,
                chart_stations,
                args.use_averages,
            )
    if args.plot_anomalies:
        chart_temperature_anomalies(inlet_list, args.show_figure)
        chart_salinity_anomalies(inlet_list, args.show_figure)
    plt.close()


if __name__ == "__main__":
    main()
