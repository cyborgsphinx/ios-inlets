import datetime
import gsw
import itertools
import logging
import numpy
import pandas
import re
from shapely.geometry import Point, Polygon
import warnings
import xarray

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

def reinsert_nan(data, placeholder):
    return numpy.fromiter(
        map(lambda x: numpy.nan if x == placeholder else x, data),
        float,
        count=len(data))

def get_array(array):
    if isinstance(array, xarray.DataArray):
        values = array.values
        return reinsert_nan(values, -99.0) if values.dtype == numpy.dtype("float64") else values
    elif isinstance(array, (float, numpy.datetime64)):
        return numpy.full(1, array)
    else:
        print("get_array called with", array)
        return array

def get_scalar(array):
    return array.item()

def find_first(source: list[str], prefix: str):
    for i, s in enumerate(source):
        if s.startswith(prefix):
            return i
    else:
        return -1

def to_float(source):
    if isinstance(source, float):
        return source
    elif isinstance(source, bytes):
        return numpy.nan if source.strip() == b"' '" else float(source.strip().decode("utf-8"))
    else:
        return numpy.nan if source.strip() == "' '" else float(source.strip())

def find_any(source, attrs: list[str]):
    for attr in attrs:
        if hasattr(source, attr):
            return getattr(source, attr)
    else:
        return None

def find_temperature_data(data):
    return find_any(data, ["TEMPRTN1", "TEMPST01", "TEMPPR01", "TEMPPR03", "TEMPS901", "TEMPS601"])

def find_salinity_data(data):
    return find_any(data, ["PSLTZZ01", "ODSDM021", "SSALST01", "PSALST01", "PSALBST1", "sea_water_practical_salinity"])

def find_oxygen_data(data):
    return find_any(data, ["DOXYZZ01", "DOXMZZ01"])

def find_depth_data(data):
    return find_any(data, ["depth", "instrument_depth", "PPSAADCP"])

def find_pressure_data(data):
    return find_any(data, ["PRESPR01", "sea_water_pressure"])

def extend_arr(arr, length):
    return numpy.full(length, arr.item()) if arr.size == 1 else arr

def is_valid_field(min_val, max_val):
    blank = min_val.strip() == "' '" or max_val.strip() == "' '"
    return blank or to_float(min_val) <= to_float(max_val)

def extract_data(source, index, min_vals, max_vals, info=None):
    if index >= 0 and is_valid_field(min_vals[index], max_vals[index]):
        replace = to_float(info["Pad"][index]) if info is not None else None
        return reinsert_nan([to_float(d[index]) for d in source], replace)
    else:
        return None

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
        pressure_dbar=None,
        filename="unknown file"):
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
                logging.warning(f"Issues arose in {filename} with {oxygen_umol_kg}, {temperature_C}, {salinity_SP}, {pressure_dbar}, assuming density to avoid invalid computations")
                assumed_density = True
                density = numpy.full(len(oxygen_umol_kg), gsw.rho([0], [0], 0)[0])
    else:
        # missing data necessary to create accurate density, assuming constant density using sigma-theta method
        logging.info(f"Not enough data in {filename} to accurately compute density. Proceeding with density as though all values are 0")
        assumed_density = True
        density = numpy.full(len(oxygen_umol_kg), gsw.rho([0], [0], 0)[0])
    return numpy.fromiter(map(lambda o, d: o * d / oxygen_umol_per_ml * metre_cube_per_litre, oxygen_umol_kg, density), float, count=len(oxygen_umol_kg)), assumed_density

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
        self.filename = filename.lower()
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

    def has_data_from(self, file_name, before=None):
        temperature_data = filter(lambda x: x.time.year < before.year, self.temperature_data) if before is not None else self.temperature_data
        salinity_data = filter(lambda x: x.time.year < before.year, self.salinity_data) if before is not None else self.salinity_data
        oxygen_data = filter(lambda x: x.time.year < before.year, self.oxygen_data) if before is not None else self.oxygen_data
        for datum in itertools.chain(temperature_data, salinity_data, oxygen_data):
            if datum.filename == file_name:
                return True
        return False

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
        longitude, latitude = None, None
        if isinstance(data, dict):
            if "longitude" in data:
                longitude = data["longitude"]
            elif "LONGITUDE" in data:
                longitude = data["LONGITUDE"]
            else:
                logging.warning("data does not contain longitude information")
                return False
            if "latitude" in data:
                latitude = data["latitude"]
            elif "LATITUDE" in data:
                latitude = data["LATITUDE"]
            else:
                logging.warning("data does not contain latitude information")
                return False
        else:
            if not hasattr(data, "longitude") or not hasattr(data, "latitude"):
                logging.warning("data does not contain full coordinate information")
                return False
            longitude = data["longitude"]
            latitude = data["latitude"]
        return self.polygon.contains(Point(longitude, latitude))

    def is_shallow(self, depth):
        return is_in_bounds(depth, *self.shallow_bounds)

    def is_middle(self, depth):
        return is_in_bounds(depth, *self.middle_bounds)

    def is_deep(self, depth):
        return is_in_bounds(depth, *self.deep_bounds)

    def produce_data(
            self,
            times,
            depths,
            data,
            longitude,
            latitude,
            filename,
            computed=False,
            assumed_density=False):
        times = extend_arr(times, len(data))
        depths = extend_arr(depths, len(data))
        if len(times) != len(data) or len(depths) != len(data):
            logging.warning(f"Data from {filename} contains times, depths, and data of different lengths")

        out = []
        once = [False] * 2
        for t, d, datum in zip(times, depths, data):
            if (not isinstance(t, datetime.datetime) and numpy.isnan(t)) or numpy.isnan(d) or numpy.isnan(datum):
                continue
            # Some data, particularly salinity data, seems to be the result of performing calculations on NaN values.
            # This data is consistently showing up as 9.96921e+36, which may relate to the "Fill Value" in creating netCDF files.
            # In any case, it appears to be as invalid as NaN, so it's being filtered out accordingly
            if datum > 9.9e+36:
                if not once[0]:
                    logging.warning(f"Data from {filename} is larger than 9.9e+36, it may have been calulated poorly")
                    once[0] = True
                continue
            if datum == -99.0:
                if not once[1]:
                    logging.warning(f"Data from {filename} has value -99, which is likely a standin for NaN")
                    once[1] = True
                continue
            if self.is_shallow(d):
                category = SHALLOW
            elif self.is_middle(d):
                category = MIDDLE
            elif self.is_deep(d):
                category = DEEP
            else:
                continue
            out.append(
                InletData(
                    get_datetime(t),
                    category,
                    datum,
                    longitude,
                    latitude,
                    filename,
                    computed=computed,
                    assumed_density=assumed_density))
        return out

    def add_temperature_data(
            self,
            data,
            depth,
            time,
            longitude,
            latitude,
            filename):
        if data is None or depth is None:
            return

        self.temperature_data.extend(self.produce_data(
            time,
            depth,
            data,
            longitude,
            latitude,
            filename))

    def add_salinity_data(
            self,
            data,
            units,
            depth,
            time,
            longitude,
            latitude,
            filename):
        if data is None or depth is None:
            return

        computed = False
        if units.lower() in ["ppt"]:
            data = gsw.SP_from_SK(data)
            computed = True
        elif units.lower() not in ["psu", "pss-78"]:
            warn_wrong_units("PSU", units, filename)
            return

        self.salinity_data.extend(self.produce_data(
            time,
            depth,
            data,
            longitude,
            latitude,
            filename,
            computed=computed))

    def add_oxygen_data(
            self,
            data,
            units,
            depth,
            time,
            temperature,
            salinity,
            pressure,
            longitude,
            latitude,
            filename):
        if data is None or depth is None:
            return

        computed = False
        assumed_density = False
        if units.lower() == "umol/kg":
            data, assumed_density = convert_umol_kg_to_mL_L(
                data,
                longitude,
                latitude,
                temperature,
                salinity,
                pressure,
                filename=filename)
            computed = True
        elif units.lower() != "ml/l":
            warn_wrong_units("mL/L", units, filename)
            return

        self.oxygen_data.extend(self.produce_data(
            time,
            depth,
            data,
            longitude,
            latitude,
            filename,
            computed=computed,
            assumed_density=assumed_density))

    def add_data_from_netcdf(self, data):
        time, longitude, latitude, filename = get_array(data.time), get_scalar(data.longitude), get_scalar(data.latitude), get_scalar(data.filename)

        depth = find_depth_data(data)
        if depth is None:
            warn_unknown_variable(data, "depth")

        temperature = find_temperature_data(data)
        if temperature is None:
            warn_unknown_variable(data, "temperature")

        salinity = find_salinity_data(data)
        if salinity is None:
            warn_unknown_variable(data, "salinity")

        oxygen = find_oxygen_data(data)
        if oxygen is None:
            warn_unknown_variable(data, "oxygen")

        pressure = find_pressure_data(data)
        if pressure is None:
            warn_unknown_variable(data, "pressure")

        self.add_temperature_data(
            get_array(temperature),
            get_array(depth),
            time,
            longitude,
            latitude,
            filename)

        self.add_salinity_data(
            get_array(salinity),
            None if salinity is None else salinity.units,
            get_array(depth),
            time,
            longitude,
            latitude,
            filename)

        self.add_oxygen_data(
            get_array(oxygen),
            None if oxygen is None else oxygen.units,
            get_array(depth),
            time,
            get_array(temperature),
            get_array(salinity),
            get_array(pressure),
            longitude,
            latitude,
            filename)

    def add_data_from_shell(self, data):
        names, units, min_vals, max_vals = data.channels["Name"], data.channels["Units"], data.channels["Minimum"], data.channels["Maximum"]

        time, longitude, latitude = numpy.full(len(data.data), data.start_dateobj), data.location["LONGITUDE"], data.location["LATITUDE"]

        depth_idx = find_first(names, "Depth")
        if depth_idx < 0:
            logging.warning(f"Shell data from {data.filename} lacks depth information. Skipping.")
            return
        depth_data = extract_data(data.data, depth_idx, min_vals, max_vals, data.channel_details)

        temperature_idx = find_first(names, "Temperature")
        temperature_data = extract_data(data.data, temperature_idx, min_vals, max_vals, data.channel_details)

        salinity_idx = find_first(names, "Salinity")
        salinity_data = extract_data(data.data, salinity_idx, min_vals, max_vals, data.channel_details)

        oxygen_idx = find_first(names, "Oxygen")
        oxygen_data = extract_data(data.data, oxygen_idx, min_vals, max_vals, data.channel_details)

        pressure_idx = find_first(names, "Pressure")
        pressure_data = extract_data(data.data, pressure_idx, min_vals, max_vals, data.channel_details)

        self.add_temperature_data(
            temperature_data,
            depth_data,
            time,
            longitude,
            latitude,
            data.filename)

        self.add_salinity_data(
            salinity_data,
            units[salinity_idx].strip(),
            depth_data,
            time,
            longitude,
            latitude,
            data.filename)

        self.add_oxygen_data(
            oxygen_data,
            units[oxygen_idx].strip(),
            depth_data,
            time,
            temperature_data,
            salinity_data,
            pressure_data,
            longitude,
            latitude,
            data.filename)
