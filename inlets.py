import json
import gsw
import logging
import numpy
import pandas
import re
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

def find_temperature_data(data):
    return find_any(data, ["TEMPRTN1", "TEMPST01", "TEMPPR01", "TEMPPR03", "TEMPS901", "TEMPS601"])

def find_salinity_data(data):
    return find_any(data, ["PSLTZZ01", "ODSDM021", "sea_water_practical_salinity"])

def find_oxygen_data(data):
    return find_any(data, ["DOXYZZ01", "DOXMZZ01"])

def find_depth_data(data):
    return find_any(data, ["depth", "instrument_depth", "PPSAADCP"])

def find_pressure_data(data):
    return find_any(data, ["PRESPR01", "sea_water_pressure"])

def warn_unknown_variable(data, var):
    # check if there is a potential variable based on the broader name
    var_list = []
    for key in data.keys():
        if re.search(var, getattr(data[key], "long_name", "").lower()):
            var_list.append(key)
    if len(var_list) != 0:
        logging.warning(f"{get_scalar(data.filename)} has unknown {var} variable. Possible values: {var_list}")

def warn_wrong_units(expected, actual, filename):
    logging.warning(f"Cowardly refusing to perform the conversion from {actual} to {expected} in {filename}")

def convert_umol_kg_to_mL_L(oxygen_umol_kg, temperature_C, salinity_SP, pressure_dbar, longitude, latitude):
    salinity_SA = gsw.SA_from_SP(salinity_SP, pressure_dbar, longitude, latitude)
    # density in kg/m^3
    density = gsw.rho(
        salinity_SA,
        gsw.CT_from_t(salinity_SA, temperature_C, pressure_dbar),
        pressure_dbar
    )
    # oxygen in umol/kg
    # conversion rate in (roughly) umol/mL
    oxygen_umol_per_ml = 44.661
    # 1 L = 10^-3 m^3
    metre_cube_per_litre = 0.001
    return list(map(lambda o, d: o * d / oxygen_umol_per_ml * metre_cube_per_litre, oxygen_umol_kg, density))

class Inlet(object):
    def __init__(self, name: str, polygon: Polygon, boundaries: list[int]):
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

    def has_temperature_data(self):
        return len(self.temperatures["shallow"]) > 0 or len(self.temperatures["middle"]) > 0 or len(self.temperatures["deep"]) > 0

    def has_salinity_data(self):
        return len(self.salinities["shallow"]) > 0 or len(self.salinities["middle"]) > 0 or len(self.salinities["deep"]) > 0

    def has_oxygen_data(self):
        return len(self.oxygens["shallow"]) > 0 or len(self.oxygens["middle"]) > 0 or len(self.oxygens["deep"]) > 0

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
        if len(times) != len(data) or len(depths) != len(data):
            logging.warning("times, depths, and data are of different lengths")
        used = False
        for t, d, datum in zip(times, depths, data):
            # Some data, particularly salinity data, seems to be the result of performing calculations on NaN values.
            # This data is consistently showing up as 9.96921e+36, which may relate to the "Fill Value" in creating netCDF files.
            # In any case, it appears to be as invalid as NaN, so it's being filtered out accordingly
            if numpy.isnan(t) or numpy.isnan(d) or numpy.isnan(datum) or datum > 9.9e+36:
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

    def add_data_to_col(self, col, time, depth, data):
        if time.size == 1:
            time = [time.item() for _ in range(len(data))]
        else:
            time = get_array(time)
        if depth.size == 1:
            depth = [depth.item() for _ in range(len(data))]
        else:
            depth = get_array(depth)
        return self.add_data(col, time, depth, data)

    def add_temperature_data_from(self, data):
        if (temp := find_temperature_data(data)) is not None:
            if (depth := find_depth_data(data)) is not None:
                return self.add_data_to_col(self.temperatures, data.time, depth, temp)
            else:
                warn_unknown_variable(data, "depth")
        else:
            warn_unknown_variable(data, "temperature")
        return False

    def add_salinity_data_from(self, data):
        if (sal := find_salinity_data(data)) is not None:
            if sal.units.lower() in ["ppt"]:
                sal = gsw.SP_from_SK(sal)
            elif sal.units.lower() not in ["psu", "pss-78"]:
                warn_wrong_units("PSU", sal.units, get_scalar(data.filename))
                return False
            if (depth := find_depth_data(data)) is not None:
                return self.add_data_to_col(self.salinities, data.time, depth, sal)
            else:
                warn_unknown_variable(data, "depth")
        else:
            warn_unknown_variable(data, "salinity")
        return False

    def add_oxygen_data_from(self, data):
        if (oxy := find_oxygen_data(data)) is not None:
            if oxy.units.lower() != "ml/l":
                if oxy.units.lower() == "umol/kg" and\
                        (temps := find_temperature_data(data)) is not None and\
                        (sal := find_salinity_data(data)) is not None and\
                        (pres := find_pressure_data(data)) is not None:
                    oxy = convert_umol_kg_to_mL_L(oxy, temps, sal, pres, get_scalar(data.longitude), get_scalar(data.latitude))
                else:
                    warn_wrong_units("mL/L", oxy.units, get_scalar(data.filename))
                    return False
            if (depth := find_depth_data(data)) is not None:
                return self.add_data_to_col(self.oxygens, data.time, depth, oxy)
            else:
                warn_unknown_variable(data, "depth")
        else:
            warn_unknown_variable(data, "oxygen")
        return False

    def add_station_from(self, data):
        year = get_datetime(data.time.head(1).item()).year
        if year not in self.stations:
            self.stations[year] = set()
        self.stations[year].add(get_scalar(data.filename))

    def add_data_from(self, data):
        if self.add_temperature_data_from(data):
            self.add_station_from(data)
        if self.add_salinity_data_from(data):
            self.add_station_from(data)
        if self.add_oxygen_data_from(data):
            self.add_station_from(data)
