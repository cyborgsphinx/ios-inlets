import gsw
import itertools
import logging
import numpy
import pandas
import re
from shapely.geometry import Point, Polygon
import warnings

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

def convert_umol_kg_to_mL_L(
        oxygen_umol_kg,
        longitude,
        latitude,
        temperature_C=None,
        salinity_SP=None,
        pressure_dbar=None):
    oxygen_umol_per_ml = 44.661
    metre_cube_per_litre = 0.001
    if temperature_C is not None and salinity_SP is not None and pressure_dbar is not None:
        assumed_density = False
        with warnings.catch_warnings(record=True) as warn:
            salinity_SA = gsw.SA_from_SP(salinity_SP, pressure_dbar, longitude, latitude)
            density = gsw.rho(
                salinity_SA,
                gsw.CT_from_t(salinity_SA, temperature_C, pressure_dbar),
                pressure_dbar
            )
            if len(warn) > 0:
                logging.warning(f"Issues arose with {oxygen_umol_kg}, {temperature_C}, {salinity_SP}, {pressure_dbar}, assuming density to avoid invalid computations")
                assumed_density = True
                density = numpy.full(len(oxygen_umol_kg), gsw.rho([0], [0], 0)[0])
    else:
        # missing data necessary to create accurate density, assuming constant density using sigma-theta method
        logging.info(f"Not enough data to accurately compute density. Proceeding with density as though all values are 0")
        assumed_density = True
        density = numpy.full(len(oxygen_umol_kg), gsw.rho([0], [0], 0)[0])
    return [o * d / oxygen_umol_per_ml * metre_cube_per_litre for o, d in zip(oxygen_umol_kg, density)], assumed_density

def get_data(col, bucket, before=None):
    data = [[datum.time, datum.datum] for datum in col if datum.bucket == bucket]
    if before is not None:
        data = [[t, d] for t, d in data if t < before]
    return zip(*data) if len(data) > 0 else [[], []]

class InletData(object):
    def __init__(
            self,
            time,
            bucket,
            datum,
            longitude,
            latitude,
            filename,
            computed=False,
            assumed_density=False):
        self.time = time
        self.bucket = bucket
        self.datum = datum
        self.longitude = longitude
        self.latitude = latitude
        self.filename = filename
        self.computed = computed
        self.assumed_density = assumed_density

class Inlet(object):
    def __init__(self, name: str, polygon: Polygon, boundaries: list[int]):
        self.name = name
        self.shallow_bounds = (boundaries[0], boundaries[1])
        self.middle_bounds = (boundaries[1], boundaries[2])
        self.deep_bounds = (boundaries[2], None)
        self.temperature_data = []
        self.salinity_data = []
        self.oxygen_data = []
        self.polygon = polygon

    def get_temperature_data(self, bucket, before=None):
        return get_data(self.temperature_data, bucket, before)

    def get_salinity_data(self, bucket, before=None):
        return get_data(self.salinity_data, bucket, before)

    def get_oxygen_data(self, bucket, before=None):
        return get_data(self.oxygen_data, bucket, before)

    def has_temperature_data(self):
        return len(self.temperature_data) > 0

    def has_salinity_data(self):
        return len(self.salinity_data) > 0

    def has_oxygen_data(self):
        return len(self.oxygen_data) > 0

    def get_station_data(self, before=None):
        temperature_data = filter(lambda x: x.time.year < before.year, self.temperature_data) if before is not None else self.temperature_data
        salinity_data = filter(lambda x: x.time.year < before.year, self.salinity_data) if before is not None else self.salinity_data
        oxygen_data = filter(lambda x: x.time.year < before.year, self.oxygen_data) if before is not None else self.oxygen_data
        stations = {}
        for datum in itertools.chain(temperature_data, salinity_data, oxygen_data):
            year = datum.time.year
            if year not in stations:
                stations[year] = set()
            stations[year].add(datum.filename)
        return stations

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

    def add_data(
            self,
            col,
            times,
            depths,
            data,
            longitude,
            latitude,
            filename,
            computed=False,
            assumed_density=False):
        if times.size == 1:
            times = numpy.full(len(data), times.item())
        else:
            times = get_array(times)
        if depths.size == 1:
            depths = numpy.full(len(data), depths.item())
        else:
            depths = get_array(depths)
        if len(times) != len(data) or len(depths) != len(data):
            logging.warning(f"Data from {filename} contains times, depths, and data of different lengths")

        once = [False] * 3
        for t, d, datum in zip(times, depths, data):
            if numpy.isnan(t) or numpy.isnan(d) or numpy.isnan(datum):
                continue
            # Some data, particularly salinity data, seems to be the result of performing calculations on NaN values.
            # This data is consistently showing up as 9.96921e+36, which may relate to the "Fill Value" in creating netCDF files.
            # In any case, it appears to be as invalid as NaN, so it's being filtered out accordingly
            if datum > 9.9e+36:
                if not once[0]:
                    logging.warning(f"Data from {filename} may have been calulated poorly")
                    once[0] = True
                continue
            if d < 0:
                if not once[1]:
                    logging.warning(f"Data from {filename} includes negative depth, and may have other incorrect data")
                    once[1] = True
                continue
            if datum == -99.0:
                if not once[2]:
                    logging.warning(f"Data from {filename} has value -99, which is likely a standin for NaN")
                    once[2] = True
                continue
            if self.is_shallow(d):
                category = SHALLOW
            elif self.is_middle(d):
                category = MIDDLE
            elif self.is_deep(d):
                category = DEEP
            else:
                continue
            col.append(
                InletData(
                    get_datetime(t),
                    category,
                    datum,
                    longitude,
                    latitude,
                    filename,
                    computed=computed,
                    assumed_density=assumed_density))

    def add_temperature_data_from(self, data):
        temperature = find_temperature_data(data)
        if temperature is None:
            warn_unknown_variable(data, "temperature")
            return
        depth = find_depth_data(data)
        if depth is None:
            warn_unknown_variable(data, "depth")
            return

        self.add_data(
            self.temperature_data,
            data.time,
            depth,
            temperature,
            data.longitude,
            data.latitude,
            get_scalar(data.filename))

    def add_salinity_data_from(self, data):
        salinity = find_salinity_data(data)
        if salinity is None:
            warn_unknown_variable(data, "salinity")
            return
        depth = find_depth_data(data)
        if depth is None:
            warn_unknown_variable(data, "depth")
            return

        computed = False
        if salinity.units.lower() in ["ppt"]:
            salinity = gsw.SP_from_SK(salinity)
            computed = True
        elif salinity.units.lower() not in ["psu", "pss-78"]:
            warn_wrong_units("PSU", salinity.units, get_scalar(data.filename))
            return

        self.add_data(
            self.salinity_data,
            data.time,
            depth,
            salinity,
            data.longitude,
            data.latitude,
            get_scalar(data.filename),
            computed=computed)

    def add_oxygen_data_from(self, data):
        oxygen = find_oxygen_data(data)
        if oxygen is None:
            warn_unknown_variable(data, "oxygen")
            return
        depth = find_depth_data(data)
        if depth is None:
            warn_unknown_variable(data, "depth")
            return

        computed = False
        assumed_density = False
        if oxygen.units.lower() != "ml/l":
            if oxygen.units.lower() == "umol/kg":
                oxygen, assumed_density = convert_umol_kg_to_mL_L(
                    oxygen,
                    get_scalar(data.longitude),
                    get_scalar(data.latitude),
                    find_temperature_data(data),
                    find_salinity_data(data),
                    find_pressure_data(data))
            else:
                warn_wrong_units("mL/L", oxygen.units, get_scalar(data.filename))
                return

        self.add_data(
            self.oxygen_data,
            data.time,
            depth,
            oxygen,
            data.longitude,
            data.latitude,
            get_scalar(data.filename),
            computed=computed,
            assumed_density=assumed_density)

    def add_data_from(self, data):
        self.add_temperature_data_from(data)
        self.add_salinity_data_from(data)
        self.add_oxygen_data_from(data)
