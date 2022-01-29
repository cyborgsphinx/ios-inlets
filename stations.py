import fnmatch
import json
import logging
import os
from shapely.geometry import Polygon

import inlets

import ios_shell.shell as ios


def main():
    shell_exts = ["bot", "che", "ctd", "ubc", "med", "xbt", "adcp", "cur"]
    # make a list of all elements in shell_exts followed by their str.upper() versions
    exts = [
        item
        for sublist in [[ext, ext.upper()] for ext in shell_exts]
        for item in sublist
    ]
    inlet_list = []
    stations = {}
    with open("inlets.geojson") as f:
        contents = json.load(f)["features"]
        for content in contents:
            name = content["properties"]["name"]
            boundaries = content["properties"]["boundaries"]
            limits = (
                content["properties"]["limits"]
                if "limits" in content["properties"]
                else {}
            )
            polygon = Polygon(content["geometry"]["coordinates"][0])
            inlet_list.append(inlets.Inlet(name, polygon, boundaries, limits))
            stations[name] = []

    for root, dirs, files in os.walk("data"):
        for ext in exts:
            for item in fnmatch.filter(files, f"*.{ext}"):
                file_name = os.path.join(root, item)
                try:
                    shell = ios.ShellFile.fromfile(file_name, process_data=False)
                except Exception as e:
                    logging.exception(f"Failed to read {file_name}: {e}")
                    continue
                for inlet in inlet_list:
                    if inlet.contains(shell.get_location()):
                        stations[inlet.name].append(
                            {
                                "type": "Feature",
                                "properties": {"station": shell.location.station},
                                "geometry": {
                                    "type": "Point",
                                    "coordinates": [
                                        shell.location.longitude,
                                        shell.location.latitude,
                                    ],
                                },
                            }
                        )
        if "HISTORY" in dirs:
            dirs.remove("HISTORY")

    for name, data in stations.items():
        print("Stations along", name)
        geo = {"type": "FeatureCollection", "features": data}
        print(json.dumps(geo))


if __name__ == "__main__":
    main()
