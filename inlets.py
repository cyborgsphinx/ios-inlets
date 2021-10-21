import json
import logging
import numpy
import pandas
from shapely.geometry import Point, Polygon

def find_polygon_for(name):
    file = "inlets.json"
    with open(file) as f:
        inlets = json.load(f)
    for inlet in inlets["features"]:
        if inlet["properties"]["name"] == name:
            return Polygon(inlet["geometry"]["coordinates"][0])
    else:
        raise RuntimeError(f"Could not find {name} in {file}")

def is_in_bounds(val, lower, upper):
    if upper is not None:
        return lower < val <= upper
    else:
        return lower < val

def get_datetime(d):
    return pandas.Timestamp(d).to_pydatetime()

def get_array(array):
    return array.values

def get_scalar(array):
    return array.item()

def find_any(source, attrs: list[str]):
    for attr in attrs:
        if hasattr(source, attr):
            return getattr(source, attr)
    else:
        return None

def warn_unknown_variable(data, var, pred):
    keys = data.keys()
    if len([key for key in keys if pred(key)]) != 0:
        logging.warning(f"{get_scalar(data.name)} has unknown {var} variable")

class Inlet(object):
    def __init__(self, name: str, polygon: Polygon, boundaries: list[int]):
        # TODO: inject polygon instead of searching for it here
        self.name = name
        self.shallow_bounds = (boundaries[0], boundaries[1])
        self.middle_bounds = (boundaries[1], boundaries[2])
        self.deep_bounds = (boundaries[2], None)
        self.temperatures = {
            "shallow": [],
            "middle": [],
            "deep": [],
        }
        self.salinities = {
            "shallow": [],
            "middle": [],
            "deep": [],
        }
        self.oxygens = {
            "shallow": [],
            "middle": [],
            "deep": []
        }
        self.stations = {}
        self.polygon = polygon

    def contains(self, data):
        longitude = data["longitude"]
        latitude = data["latitude"]
        return self.polygon.contains(Point(longitude, latitude))

    def is_shallow(self, depth):
        return is_in_bounds(depth, *self.shallow_bounds)

    def is_middle(self, depth):
        return is_in_bounds(depth, *self.middle_bounds)

    def is_deep(self, depth):
        return is_in_bounds(depth, *self.deep_bounds)

    def add_data(self, col, times, depths, data):
        if len(times) != len(depths):
            logging.warning("times and depths are of different lengths")
        used = False
        for t, d, datum in zip(times, depths, data):
            if numpy.isnan(t) or numpy.isnan(d) or numpy.isnan(datum):
                continue
            if self.is_shallow(d):
                category = "shallow"
            elif self.is_middle(d):
                category = "middle"
            elif self.is_deep(d):
                category = "deep"
            else:
                continue
            col[category].append([get_datetime(t), datum])
            used = True
        return used

    def add_data_constant_time(self, col, time, depths, data):
        return self.add_data(col, [time for _ in range(len(depths))], depths, data)

    def add_data_constant_depth(self, col, times, depth, data):
        return self.add_data(col, times, [depth for _ in range(len(times))], data)

    def add_temperatures(self, times, depths, temperatures):
        return self.add_data(self.temperatures, times, depths, temperatures)

    def add_temperatures_constant_time(self, time, depths, temperatures):
        return self.add_data_constant_time(self.temperatures, time, depths, temperatures)

    def add_temperatures_constant_depth(self, times, depth, temperatures):
        return self.add_data_constant_depth(self.temperatures, times, depth, temperatures)

    def add_salinities(self, times, depths, salinities):
        return self.add_data(self.salinities, times, depths, salinities)

    def add_salinities_constant_time(self, time, depths, salinities):
        return self.add_data_constant_time(self.salinities, time, depths, salinities)

    def add_salinities_constant_depth(self, times, depth, salinities):
        return self.add_data_constant_depth(self.salinities, times, depth, salinities)

    def add_oxygen_data(self, times, depths, oxygen_data):
        return self.add_data(self.oxygens, times, depths, oxygen_data)

    def add_oxygen_data_constant_time(self, time, depths, oxygen_data):
        return self.add_data_constant_time(self.oxygens, time, depths, oxygen_data)

    def add_oxygen_data_constant_depth(self, times, depth, oxygen_data):
        return self.add_data_constant_time(self.oxygens, times, depth, oxygen_data)

    def add_temperature_data_from(self, data):
        if not hasattr(data, "depth"):
            logging.info(f"data from {data.filename.item()} has no depth information, discarding")
            return False
        # find values in specific depth intervals
        # will be plotted against time
        # BOT/1930-031-0001.bot.nc and CTD/1966-062-0129.ctd.nc used as example
        if (temps := find_any(data, ["TEMPRTN1", "TEMPST01"])) is not None:
            depth = get_array(data.depth)
            if data.time.size == 1:
                time = get_scalar(data.time)
                return self.add_temperatures_constant_time(time, depth, temps)
            else:
                time = get_array(data.time)
                return self.add_temperatures(time, depth, temps)
        # ADCP/nep1_20060512_20060525_0095m.adcp.L1.nc and CUR/CM1_19890407_19890504_0020m.cur.nc used as example
        elif (temps := find_any(data, ["TEMPPR01", "TEMPPR03"])) is not None:
            time = get_array(data.time)
            if hasattr(data, "PPSAADCP"):
                # treat like ADCP/nep1_20060512_20060525_0095m.adcp.L1.nc
                depth = get_array(data.PPSAADCP)
                return self.add_temperatures(time, depth, temps)
            elif hasattr(data, "instrument_depth"):
                # treat like CUR/CM1_19890407_19890504_0020m.cur.nc
                return self.add_temperatures_constant_depth(time, data.instrument_depth.item(), temps)
            else:
                name = getattr(data, "filename", None)
                file_name = get_scalar(name) if name is not None else "unknown file"
                logging.warning(f"{file_name} has unknown depth variable")
        # CTD/2021-020-0001.ctd.nc and BOT/1983-030-0018.che.nc used as example
        elif (temps := find_any(data, ["TEMPS901", "TEMPS601"])) is not None:
            depth = get_array(data.depth)
            time = get_scalar(data.time)
            return self.add_temperatures_constant_time(time, depth, temps)
        else:
            warn_unknown_variable(data, "temperature", lambda s: s.lower().startswith("temp"))
        return False

    def add_salinity_data_from(self, data):
        # TODO: finish
        # CTD/1966-062-0129.ctd.nc used as example
        #if (sal := find_any(data, ["SSALST01"])) is not None:
        #    depth = get_array(data.depth)
        #    time = get_scalar(data.time)
        #    return self.add_salinities_constant_time(time, depth, sal)
        #else:
        #    warn_unknown_variable(data, "salinity", lambda s: s.lower().startswith("ssal"))
        return False

    def add_oxygen_data_from(self, data):
        # TODO: finish
        # BOT/1978-033-0013.bot.nc used as example
        #if (oxy := find_any(data, ["DOXYZZ01"])) is not None:
        #    depth = get_array(data.depth)
        #    time = get_scalar(data.time)
        #    return self.add_oxygen_data_constant_time(time, depth, oxy)
        #else:
        #    warn_unknown_variable(data, "oxygen", lambda s: s.lower().startswith("doxy"))
        return False

    def add_station_from(self, data):
        year = get_datetime(data.time.head(1).item()).year
        if year not in self.stations:
            self.stations[year] = set()
        self.stations[year].add(get_scalar(data.filename))

    def add_data_from(self, data):
        if self.add_temperature_data_from(data):
            self.add_station_from(data)
