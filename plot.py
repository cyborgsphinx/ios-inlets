import argparse
import datetime
import inlets
import itertools
import matplotlib
import matplotlib.pyplot as plt
import os
from typing import Dict, List

END = datetime.datetime.now()
INLET_LINES = ["m-s", "y-d", "k-o", "c-^", "b-d", "g-s", "r-s", "m-d"]


###################
# Utility functions
###################


def normalize(string: str):
    return string.strip().lower().replace(" ", "-")


def figure_path(filename: str):
    return os.path.join("figures", filename)


def label_from_bounds(lower, upper):
    if upper is None:
        return f">{lower}m"
    else:
        return f"{lower}m-{upper}m"


########################
# Single inlet functions
########################


def chart_deep_data(inlet: inlets.Inlet, limits: List[float], data_fn):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    shallow_time, shallow_data = data_fn(inlet, inlets.DEEP)
    middle_time, middle_data = data_fn(inlet, inlets.DEEPER)
    deep_time, deep_data = data_fn(inlet, inlets.DEEPEST)
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
        label=label_from_bounds(*inlet.deep_bounds),
    )
    plt.plot(
        middle_time,
        middle_data,
        "+m",
        label=label_from_bounds(*inlet.deeper_bounds),
    )
    plt.plot(deep_time, deep_data, "xb", label=label_from_bounds(*inlet.deepest_bounds))
    plt.legend()


def chart_surface_data(inlet: inlets.Inlet, limits: List[float], data_fn):
    plt.clf()
    surface_time, surface_data = data_fn(inlet, inlets.SURFACE)
    shallow_time, shallow_data = data_fn(inlet, inlets.SHALLOW)
    if len(limits) > 1:
        surface_time, surface_data = zip(
            *[
                [t, d]
                for t, d in zip(surface_time, surface_data)
                if limits[0] < d < limits[1]
            ]
        )
        shallow_time, shallow_data = zip(
            *[
                [t, d]
                for t, d in zip(shallow_time, shallow_data)
                if limits[0] < d < limits[1]
            ]
        )
    plt.plot(surface_time, surface_data, "xg", label=label_from_bounds(*inlet.surface_bounds))
    if inlet.shallow_bounds is not None:
        plt.plot(shallow_time, shallow_data, "+m", label=label_from_bounds(*inlet.shallow_bounds))
    plt.legend()


def chart_temperatures(inlet: inlets.Inlet, limits: Dict[str, List[float]], use_averages: bool):
    average = "-average" if use_averages else ""
    ylabel = "Temperature (C)"
    data_fn = lambda inlet, bucket: inlet.get_temperature_data(
        bucket, before=END, do_average=use_averages
    )

    chart_deep_data(inlet, limits["deep"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Deep Water Temperature")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-deep-temperature{average}.png"))

    chart_surface_data(inlet, limits["surface"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Surface Water Temperature")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-surface-temperature{average}.png"))


def chart_salinities(inlet: inlets.Inlet, limits: Dict[str,List[float]], use_averages: bool):
    average = "-average" if use_averages else ""
    ylabel = "Salinity (PSU)"
    data_fn = lambda inlet, bucket: inlet.get_salinity_data(
        bucket, before=END, do_average=use_averages
    )

    chart_deep_data(inlet, limits["deep"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Deep Water Salinity")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-deep-salinity{average}.png"))

    chart_surface_data(inlet, limits["surface"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Surface Water Salinity")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-surface-salinity{average}.png"))


def chart_oxygen_data(inlet: inlets.Inlet, limits: Dict[str,List[float]], use_averages: bool):
    average = "-average" if use_averages else ""
    ylabel = "DO (ml/l)"
    data_fn = lambda inlet, bucket: inlet.get_oxygen_data(
        bucket, before=END, do_average=use_averages
    )

    chart_deep_data(inlet, limits["deep"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Deep Water Dissolved Oxygen")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-deep-oxygen{average}.png"))

    chart_surface_data(inlet, limits["surface"], data_fn)
    plt.ylabel(ylabel)
    plt.title(f"{inlet.name} Surface Water Dissolved Oxygen")
    plt.savefig(figure_path(f"{normalize(inlet.name)}-surface-oxygen{average}.png"))


def chart_stations(inlet: inlets.Inlet, _limits: Dict[str,List[float]], _use_averages: bool):
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
    plt.savefig(figure_path(f"{normalize(inlet.name)}-samples.png"))


def do_chart(
    inlet: inlets.Inlet,
    kind: str,
    use_limits: bool,
    chart_fn,
    use_averages: bool,
):
    averaging = " averages" if use_averages else ""
    print(f"Producing {kind}{averaging} plot for {inlet.name}")
    chart_fn(
        inlet,
        inlet.limits[kind] if use_limits and kind in inlet.limits else {"deep": [], "surface": []},
        use_averages,
    )


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


def do_annual_work(inlet_list, data_fn, averaging_fn, limit_fn):
    plt.clf()
    for inlet, line_style in zip(inlet_list, INLET_LINES):
        limits = limit_fn(inlet)
        totals = {}
        times, data = data_fn(inlet)
        for time, datum in zip(times, data):
            if len(limits) > 0 and not (limits[0] < datum < limits[1]):
                continue
            year = time.year
            if year not in totals:
                totals[year] = (0, 0)
            total, num = totals[year]
            totals[year] = (total + datum, num + 1)

        averages = averaging_fn(totals, data)
        years, values = zip(*sorted(averages.items(), key=lambda item: item[0]))
        plt.plot(years, values, line_style, label=inlet.name)

    plt.legend()


def do_annual_work_parts(inlet_list, data_fn, averaging_fn, y_label, title, limit_fn):
    areas = set(inlet.area for inlet in inlet_list)
    for area in areas:
        do_annual_work([inlet for inlet in inlet_list if inlet.area == area], data_fn, averaging_fn, limit_fn)
        plt.ylabel(y_label)
        plt.title(f"{title} - {area}")
        plt.savefig(figure_path(f"{normalize(title)}-{normalize(area)}.png"))


def annual_averaging(totals, _data):
    del _data
    return {y: t / n for y, (t, n) in totals.items()}


def anomalies_averaging(totals, data):
    avg = sum(data) / len(data)
    avgs = annual_averaging(totals, data)
    return {y: a - avg for y, a in avgs.items()}


def chart_anomalies(inlet_list: List[inlets.Inlet], data_fn, y_label, title, use_limits):
    do_annual_work_parts(inlet_list, data_fn, anomalies_averaging, y_label, title, use_limits)


def chart_temperature_anomalies(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing temperature anomaly plots")
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_temperature_data(inlets.USED_DEEP, do_average=True),
        "Temperature (C)",
        "Deep Water Temperature Anomalies",
        lambda inlet: inlet.limits["temperature"]["deep"] if use_limits and "temperature" in inlet.limits else [],
    )
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_temperature_data(inlets.USED_SURFACE, do_average=True),
        "Temperature (C)",
        "Surface Water Temperature Anomalies",
        lambda inlet: inlet.limits["temperature"]["surface"] if use_limits and "temperature" in inlet.limits else [],
    )


def chart_salinity_anomalies(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing salinity anomaly plots")
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_salinity_data(inlets.USED_DEEP, do_average=True),
        "Salinity (PSU)",
        "Deep Water Salinity Anomalies",
        lambda inlet: inlet.limits["salinity"]["deep"] if use_limits and "salinity" in inlet.limits else [],
    )
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_salinity_data(inlets.USED_SURFACE, do_average=True),
        "Salinity (PSU)",
        "Surface Water Salinity Anomalies",
        lambda inlet: inlet.limits["salinity"]["surface"] if use_limits and "salinity" in inlet.limits else [],
    )


def chart_oxygen_anomalies(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing oxygen anomaly plot")
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_oxygen_data(inlets.USED_DEEP, do_average=True),
        "Oxygen (mL/L)",
        "Deep Water Dissolved Oxygen Anomalies",
        lambda inlet: inlet.limits["oxygen"]["deep"] if use_limits and "oxygen" in inlet.limits else [],
    )
    chart_anomalies(
        inlet_list,
        lambda inlet: inlet.get_oxygen_data(inlets.USED_SURFACE, do_average=True),
        "Oxygen (mL/L)",
        "Surface Water Dissolved Oxygen Anomalies",
        lambda inlet: inlet.limits["oxygen"]["surface"] if use_limits and "oxygen" in inlet.limits else [],
    )


def do_chart_annual_averages(inlet_list: List[inlets.Inlet], data_fn, y_label, title, limit_fn):
    do_annual_work_parts(inlet_list, data_fn, annual_averaging, y_label, title, limit_fn)


def chart_annual_temperature_averages(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing annual temperature plots")
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_temperature_data(inlets.USED_DEEP, do_average=True),
        "Temperature (C)",
        "Deep Water Temperature Annual Averages",
        lambda inlet: inlet.limits["temperature"]["deep"] if use_limits and "temperature" in inlet.limits else [],
    )
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_temperature_data(inlets.USED_SURFACE, do_average=True),
        "Temperature (C)",
        "Surface Water Temperature Annual Averages",
        lambda inlet: inlet.limits["temperature"]["surface"] if use_limits and "temperature" in inlet.limits else [],
    )


def chart_annual_salinity_averages(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing annual salinity plots")
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_salinity_data(inlets.USED_DEEP, do_average=True),
        "Salinity (PSU)",
        "Deep Water Salinity Annual Averages",
        lambda inlet: inlet.limits["salinity"]["deep"] if use_limits and "salinity" in inlet.limits else [],
    )
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_salinity_data(inlets.USED_SURFACE, do_average=True),
        "Salinity (PSU)",
        "Surface Water Salinity Annual Averages",
        lambda inlet: inlet.limits["salinity"]["surface"] if use_limits and "salinity" in inlet.limits else [],
    )


def chart_annual_oxygen_averages(inlet_list: List[inlets.Inlet], use_limits: bool):
    print("Producing annual oxygen plots")
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_oxygen_data(inlets.USED_DEEP, do_average=True),
        "Oxygen (mL/L)",
        "Deep Water Dissolved Oxygen Annual Averages",
        lambda inlet: inlet.limits["oxygen"]["deep"] if use_limits and "oxygen" in inlet.limits else [],
    )
    do_chart_annual_averages(
        inlet_list,
        lambda inlet: inlet.get_oxygen_data(inlets.USED_SURFACE, do_average=True),
        "Oxygen (mL/L)",
        "Surface Water Dissolved Oxygen Annual Averages",
        lambda inlet: inlet.limits["oxygen"]["surface"] if use_limits and "oxygen" in inlet.limits else [],
    )


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
        inlet.data.get_temperature_data((None, None)),
        inlet.data.get_salinity_data((None, None)),
        inlet.data.get_oxygen_data((None, None)),
    ):
        year = datum.time.year
        min_year = min(year, min_year)
        max_year = max(year, max_year)
        month = months[datum.time.month]
        if year not in files[month]:
            files[month][year] = set()
        files[month][year].add(datum.source)
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
    matplotlib.rc("axes", titlesize=25)
    matplotlib.rc("xtick", labelsize=20)
    matplotlib.rc("ytick", labelsize=20)
    plt.figure(figsize=(40, 10), constrained_layout=True)
    plt.imshow(list(map(list, zip(*values))), vmin=0, vmax=biggest, cmap="Blues")
    plt.yticks(ticks=range(12), labels=months[1:])
    plt.xticks(
        ticks=range(0, year_range, 2),
        labels=range(min_year, max_year + 1, 2),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    for i, _ in enumerate(values):
        for j, _ in enumerate(values[i]):
            plt.text(i, j, values[i][j], ha="center", va="center", color="k", fontsize="large")

    plt.title(f"{inlet.name} Sampling Frequency by Month")
    plt.axis("tight")
    plt.colorbar()
    plt.savefig(figure_path(f"{normalize(inlet.name)}-monthly-sampling.png"))
    plt.close()

    # reset values
    matplotlib.rcdefaults()
    plt.axis("auto")


def main():
    parser = argparse.ArgumentParser()
    # inlet retrieval args
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-n", "--from-netcdf", action="store_true")
    parser.add_argument("-e", "--from-erddap", action="store_true")
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    # plot args
    parser.add_argument("-l", "--no-limits", action="store_true")
    parser.add_argument("-i", "--inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-k", "--limit-name", type=str, nargs="+", default=[])
    parser.add_argument("-I", "--remove-inlet-name", type=str, nargs="+", default=[])
    parser.add_argument("-b", "--plot-buckets", action="store_true")
    parser.add_argument("-a", "--plot-averages", action="store_true")
    parser.add_argument("-A", "--plot-annual", action="store_true")
    parser.add_argument("-s", "--plot-sampling", action="store_true")
    parser.add_argument("--plot-all", action="store_true")
    args = parser.parse_args()
    inlet_list = inlets.get_inlets(
        args.data,
        from_saved=args.from_saved,
        from_netcdf=args.from_netcdf,
        from_erddap=args.from_erddap,
        inlet_names=args.inlet_name,
        drop_names=args.remove_inlet_name,
        keep_names=args.limit_name,
    )
    plt.figure(figsize=(8, 6))
    if args.plot_all:
        (plot_annual, plot_sampling, plot_average, plot_raw, plot_buckets) = (True, True, True, True, False)
    else:
        (plot_annual, plot_sampling, plot_average, plot_raw, plot_buckets) = (args.plot_annual, args.plot_sampling, args.plot_averages, args.plot_raw, args.plot_buckets)
    if plot_annual:
        chart_annual_temperature_averages(inlet_list, not args.no_limits)
        chart_annual_salinity_averages(inlet_list, not args.no_limits)
        chart_annual_oxygen_averages(inlet_list, not args.no_limits)
        chart_temperature_anomalies(inlet_list, not args.no_limits)
        chart_salinity_anomalies(inlet_list, not args.no_limits)
        chart_oxygen_anomalies(inlet_list, not args.no_limits)
    if plot_sampling:
        for inlet in inlet_list:
            do_chart(
                inlet,
                "samples",
                False,
                chart_stations,
                False,
            )
            chart_monthly_sample(inlet)
    if plot_buckets:
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.DEEP,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.DEEPER,
            chart_all_temperature,
        )
        do_chart_all(
            inlet_list,
            "temperature",
            inlets.DEEPEST,
            chart_all_temperature,
        )
        do_chart_all(inlet_list, "salinity", inlets.DEEP, chart_all_salinity)
        do_chart_all(inlet_list, "salinity", inlets.DEEPER, chart_all_salinity)
        do_chart_all(inlet_list, "salinity", inlets.DEEPEST, chart_all_salinity)
        do_chart_all(inlet_list, "oxygen", inlets.DEEP, chart_all_oxygen)
        do_chart_all(inlet_list, "oxygen", inlets.DEEPER, chart_all_oxygen)
        do_chart_all(inlet_list, "oxygen", inlets.DEEPEST, chart_all_oxygen)
    if plot_average:
        for inlet in inlet_list:
            do_chart(
                inlet,
                "temperature",
                not args.no_limits,
                chart_temperatures,
                True,
            )
            do_chart(
                inlet,
                "salinity",
                not args.no_limits,
                chart_salinities,
                True,
            )
            do_chart(
                inlet,
                "oxygen",
                not args.no_limits,
                chart_oxygen_data,
                True,
            )
    if plot_raw:
        for inlet in inlet_list:
            do_chart(
                inlet,
                "temperature",
                not args.no_limits,
                chart_temperatures,
                False,
            )
            do_chart(
                inlet,
                "salinity",
                not args.no_limits,
                chart_salinities,
                False,
            )
            do_chart(
                inlet,
                "oxygen",
                not args.no_limits,
                chart_oxygen_data,
                False,
            )
    plt.close()


if __name__ == "__main__":
    main()
