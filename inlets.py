import json
import gsw
import logging
import numpy
import pandas
import re
from shapely.geometry import Point, Polygon

SHALLOW = "shallow"
MIDDLE = "middle"
DEEP = "deep"

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

def convert_umol_kg_to_mL_L(oxygen_umol_kg, longitude, latitude, temperature_C=None, salinity_SP=None, pressure_dbar=None):
    # conversion rate in (roughly) umol/mL
    oxygen_umol_per_ml = 44.661
    # 1 L = 10^-3 m^3
    metre_cube_per_litre = 0.001
    if temperature_C is not None and salinity_SP is not None and pressure_dbar is not None:
        salinity_SA = gsw.SA_from_SP(salinity_SP, pressure_dbar, longitude, latitude)
        # density in kg/m^3
        density = gsw.rho(
            salinity_SA,
            gsw.CT_from_t(salinity_SA, temperature_C, pressure_dbar),
            pressure_dbar
        )
        # oxygen in umol/kg
        return [o * d / oxygen_umol_per_ml * metre_cube_per_litre for o, d in zip(oxygen_umol_kg, density)], False
    else:
        # missing data necessary to create accurate density, assuming constant density using sigma-theta method
        density = gsw.rho(
            [0],
            [0],
            0
        )[0]
        return [o * density / oxygen_umol_per_ml * metre_cube_per_litre for o in oxygen_umol_kg], True

class InletData(object):
    def __init__(self, time, bucket, datum, longitude, latitude, computed=False, assumed_density=False):
        self.time = time
        self.bucket = bucket
        self.datum = datum
        self.longitude = longitude
        self.latitude = latitude
        self.computed = computed
        self.assumed_density = assumed_density

class Inlet(object):
    def __init__(self, name: str, polygon: Polygon, boundaries: list[int]):
        self.name = name
        self.shallow_bounds = (boundaries[0], boundaries[1])
        self.middle_bounds = (boundaries[1], boundaries[2])
        self.deep_bounds = (boundaries[2], None)
        self.temperatures = []
        self.salinities = []
        self.oxygens = []
        self.stations = {}
        self.polygon = polygon

    def get_temperature_data(self, bucket, before=None):
        data = [[temp.time, temp.datum] for temp in self.temperatures if temp.bucket == bucket]
        if before is not None:
            data = [[t, d] for t, d in data if t < before]
        return zip(*data)

    def get_salinity_data(self, bucket, before=None):
        data = [[temp.time, temp.datum] for temp in self.salinities if temp.bucket == bucket]
        if before is not None:
            data = [[t, d] for t, d in data if t < before]
        return zip(*data)

    def get_oxygen_data(self, bucket, before=None):
        data = [[temp.time, temp.datum] for temp in self.oxygens if temp.bucket == bucket]
        if before is not None:
            data = [[t, d] for t, d in data if t < before]
        return zip(*data)

    def has_temperature_data(self):
        return len(self.temperatures) > 0

    def has_salinity_data(self):
        return len(self.salinities) > 0

    def has_oxygen_data(self):
        return len(self.oxygens) > 0

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

    def add_data(self, col, times, depths, data, longitude, latitude, computed=False, assumed_density=False):
        if times.size == 1:
            times = [times.item() for _ in range(len(data))]
        else:
            times = get_array(times)
        if depths.size == 1:
            depths = [depths.item() for _ in range(len(data))]
        else:
            depths = get_array(depths)
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
                category = SHALLOW
            elif self.is_middle(d):
                category = MIDDLE
            elif self.is_deep(d):
                category = DEEP
            else:
                continue
            col.append(InletData(get_datetime(t), category, datum, longitude, latitude, computed=computed, assumed_density=assumed_density))
            used = True
        return used

    def add_temperature_data_from(self, data):
        if (temp := find_temperature_data(data)) is not None:
            if (depth := find_depth_data(data)) is not None:
                return self.add_data(self.temperatures, data.time, depth, temp, data.longitude, data.latitude)
            else:
                warn_unknown_variable(data, "depth")
        else:
            warn_unknown_variable(data, "temperature")
        return False

    def add_salinity_data_from(self, data):
        if (sal := find_salinity_data(data)) is not None:
            computed = False
            if sal.units.lower() in ["ppt"]:
                sal = gsw.SP_from_SK(sal)
                computed = True
            elif sal.units.lower() not in ["psu", "pss-78"]:
                warn_wrong_units("PSU", sal.units, get_scalar(data.filename))
                return False
            if (depth := find_depth_data(data)) is not None:
                return self.add_data(
                    self.salinities,
                    data.time,
                    depth,
                    sal,
                    data.longitude,
                    data.latitude,
                    computed=computed)
            else:
                warn_unknown_variable(data, "depth")
        else:
            warn_unknown_variable(data, "salinity")
        return False

    def add_oxygen_data_from(self, data):
        if (oxy := find_oxygen_data(data)) is not None:
            computed = False
            assumed_density = False
            if oxy.units.lower() != "ml/l":
                if oxy.units.lower() == "umol/kg":
                    oxy, assumed_density = convert_umol_kg_to_mL_L(
                        oxy,
                        get_scalar(data.longitude),
                        get_scalar(data.latitude),
                        find_temperature_data(data),
                        find_salinity_data(data),
                        find_pressure_data(data))
                else:
                    warn_wrong_units("mL/L", oxy.units, get_scalar(data.filename))
                    return False
            if (depth := find_depth_data(data)) is not None:
                return self.add_data(
                    self.oxygens,
                    data.time,
                    depth,
                    oxy,
                    data.longitude,
                    data.latitude,
                    computed=computed,
                    assumed_density=assumed_density)
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
