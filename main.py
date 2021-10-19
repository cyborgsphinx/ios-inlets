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

PICKLE_NAME = "saanich.pickle"

def truncate(data):
    return [[d, t] for d, t in data if d < datetime.datetime(2000, 1, 1)]

def chart_temperatures(inlet):
    # produce a matplotlib chart, which can be shown or saved at the upper level
    ax = plt.axes(label=f"{inlet.name} Deep Water Temperatures")
    shallow_x, shallow_y = zip(*truncate(inlet.temperatures["shallow"]))
    ax.plot(shallow_x, shallow_y, "xg", label=f"{inlet.shallow_bounds[0]}m-{inlet.shallow_bounds[1]}m")
    middle_x, middle_y = zip(*truncate(inlet.temperatures["middle"]))
    ax.plot(middle_x, middle_y, "+m", label=f"{inlet.middle_bounds[0]}m-{inlet.middle_bounds[1]}m")
    deep_x, deep_y = zip(*truncate(inlet.temperatures["deep"]))
    ax.plot(deep_x, deep_y, "xb", label=f">{inlet.deep_bounds[0]}m")
    ax.set_ylabel("Temperature (C)")
    ax.legend()
    ax.set_title(f"{inlet.name} Deep Water Temperatures")
    return ax

def normalize(string):
    return string.strip().lower()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-s", "--show-figure", action="store_true")
    args = parser.parse_args()
    if args.from_saved:
        with open(PICKLE_NAME, mode="rb") as f:
            saanich = pickle.load(f)
    else:
        with open("inlets2.json") as f:
            contents = json.load(f)["features"][0]
        poly = Polygon(contents["geometry"]["coordinates"][0])
        saanich = inlets.Inlet(contents["properties"]["name"], poly, contents["properties"]["boundaries"])
        for root, _, files in os.walk("data"):
            #print(root, "-", dirs)
            for item in fnmatch.filter(files, "*.nc"):
                #print(". .", item)
                file_name = os.path.join(root, item)
                data = xarray.open_dataset(file_name)
                if saanich.contains(data):
                    #print("...", item)
                    try:
                        saanich.add_temperature_data_from(data)
                    except:
                        logging.exception(f"Exception occurred in {file_name}")
                        raise
                elif normalize(data["geographic_area"].item()) in ["saanich-inlet", "satellite-channel", "cowichan-bay"]:
                    logging.warning(f"{file_name} claims to be in saanich inlet, but its location is outside the polygon")
        with open(PICKLE_NAME, mode="wb") as f:
            pickle.dump(saanich, f)
    chart_temperatures(saanich)
    if args.show_figure:
        plt.show()
    else:
        plt.savefig("saanich-inlet-temperature.png")

if __name__ == "__main__":
    main()
