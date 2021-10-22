import argparse
import datetime
import inlets
import json
import fnmatch
import logging
import matplotlib.pyplot as plt
import pickle
import os
from shapely.geometry import Polygon
import xarray

PICKLE_NAME = "inlets.pickle"

def truncate(data):
    return [[d, t] for d, t in data if d < datetime.datetime(2000, 1, 1)]

def chart_temperatures(inlet):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    shallow_x, shallow_y = zip(*truncate(inlet.temperatures["shallow"]))
    plt.plot(shallow_x, shallow_y, "xg", label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m")
    middle_x, middle_y = zip(*truncate(inlet.temperatures["middle"]))
    plt.plot(middle_x, middle_y, "+m", label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m")
    deep_x, deep_y = zip(*truncate(inlet.temperatures["deep"]))
    plt.plot(deep_x, deep_y, "xb", label=f">{inlet.deep_bounds[0]}m")
    plt.ylabel("Temperature (C)")
    plt.legend()
    plt.title(f"{inlet.name} Deep Water Temperatures")

def chart_salinities(inlet):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    shallow_x, shallow_y = zip(*truncate(inlet.salinities["shallow"]))
    plt.plot(shallow_x, shallow_y, "xg", label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m")
    middle_x, middle_y = zip(*truncate(inlet.salinities["middle"]))
    plt.plot(middle_x, middle_y, "+m", label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m")
    deep_x, deep_y = zip(*truncate(inlet.salinities["deep"]))
    plt.plot(deep_x, deep_y, "xb", label=f">{inlet.deep_bounds[0]}m")
    plt.ylabel("Salinity (PSU)")
    plt.legend()
    plt.title(f"{inlet.name} Deep Water Salinity")

def chart_oxygen_data(inlet):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    plt.clf()
    shallow_x, shallow_y = zip(*truncate(inlet.oxygens["shallow"]))
    plt.plot(shallow_x, shallow_y, "xg", label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m")
    middle_x, middle_y = zip(*truncate(inlet.oxygens["middle"]))
    plt.plot(middle_x, middle_y, "+m", label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m")
    deep_x, deep_y = zip(*truncate(inlet.oxygens["deep"]))
    plt.plot(deep_x, deep_y, "xb", label=f">{inlet.deep_bounds[0]}m")
    plt.ylabel("DO (ml/l)")
    plt.legend()
    plt.title(f"{inlet.name} Deep Water Dissolved Oxygen")

def chart_stations(inlet):
    plt.clf()
    data = []
    for year, stations in inlet.stations.items():
        if year > 2000:
            continue
        data.extend([year for _ in range(len(stations))])
    n_bins = max(data) - min(data) + 1
    plt.hist(data, bins=n_bins, align="left", label=f"Number of files {len(data)}")
    plt.ylabel("Number of Stations")
    plt.legend()
    plt.title(f"{inlet.name} Sampling History")

def normalize(string):
    return string.strip().lower().replace(' ', '-')

def do_chart(inlet, kind: str, show_figure: bool, chart_fn):
    chart_fn(inlet)
    if show_figure:
        plt.show()
    else:
        plt.savefig(os.path.join("figures", f"{normalize(inlet.name)}-{kind}.png"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-s", "--show-figure", action="store_true")
    args = parser.parse_args()
    inlet_s = []
    if args.from_saved:
        with open(PICKLE_NAME, mode="rb") as f:
            inlet_s = pickle.load(f)
    else:
        with open("inlets.json") as f:
            contents = json.load(f)["features"]
            for content in contents:
                name = content["properties"]["name"]
                boundaries = content["properties"]["boundaries"]
                polygon = Polygon(content["geometry"]["coordinates"][0])
                inlet_s.append(inlets.Inlet(name, polygon, boundaries))
        for root, _, files in os.walk("data"):
            #print(root, "-", dirs)
            for item in fnmatch.filter(files, "*.nc"):
                #print(". .", item)
                file_name = os.path.join(root, item)
                data = xarray.open_dataset(file_name)
                for inlet in inlet_s:
                    if inlet.contains(data):
                        #print("...", item)
                        try:
                            inlet.add_data_from(data)
                        except:
                            logging.exception(f"Exception occurred in {file_name}")
                            raise
        with open(PICKLE_NAME, mode="wb") as f:
            pickle.dump(inlet_s, f)
    for inlet in inlet_s:
        do_chart(inlet, "temperature", args.show_figure, chart_temperatures)
        do_chart(inlet, "salinity", args.show_figure, chart_salinities)
        do_chart(inlet, "oxygen", args.show_figure, chart_oxygen_data)
        do_chart(inlet, "stations", args.show_figure, chart_stations)

if __name__ == "__main__":
    main()
