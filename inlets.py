import datetime
import gsw
import itertools
import logging
import numpy
import os
import pandas
import re
from shapely.geometry import Point, Polygon
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

def reinsert_nan(data, placeholder, length=None):
    return numpy.fromiter(
        (numpy.nan if x == placeholder else x for x in data),
        float,
        count=len(data) if length is None else length)

def get_scalar(array):
    return array.item()

def get_array(array):
    if isinstance(array, xarray.DataArray):
        values = array.values
        if values.size == 1:
            values = numpy.full(1, get_scalar(values))
        return reinsert_nan(values, -99.0) if values.dtype.kind == "f" else values
    elif isinstance(array, (float, numpy.datetime64)):
        return numpy.full(1, array)
    else:
        return array

def find_first(source: list[str], *args):
    for i, s in enumerate(source):
        for prefix in args:
            if s.startswith(prefix):
                return i
    else:
        return -1

def to_float(source):
    if isinstance(source, float):
        return source
    elif isinstance(source, bytes):
        return numpy.nan if source.strip() in [b"' '", b"n/a", b""] else float(source.strip().decode("utf-8"))
    else:
        return numpy.nan if source.strip() in ["' '", "n/a", ""] else float(source.strip())

def find_any(source, *attrs):
    for attr in attrs:
        if hasattr(source, attr):
            return getattr(source, attr)
    else:
        return None

def find_temperature_data(data):
    return find_any(data, "TEMPRTN1", "TEMPST01", "TEMPPR01", "TEMPPR03", "TEMPS901", "TEMPS601")

def find_salinity_data(data):
    return find_any(data, "PSLTZZ01", "ODSDM021", "SSALST01", "PSALST01", "PSALBST1", "sea_water_practical_salinity")

def find_oxygen_data(data):
    return find_any(data, "DOXYZZ01", "DOXMZZ01")

def find_depth_data(data):
    return find_any(data, "depth", "instrument_depth", "PPSAADCP")

def find_pressure_data(data):
    return find_any(data, "PRESPR01", "sea_water_pressure")

def extend_arr(arr, length):
    return numpy.full(length, arr.item()) if arr.size == 1 else arr

def is_valid_field(min_val, max_val):
    blank = min_val.strip() == "' '" or max_val.strip() == "' '"
    return blank or to_float(min_val) <= to_float(max_val)

def get_pad_value(info, index):
    if index < 0 or info is None:
        return None
    return to_float(info["Pad"][index])

def has_quality(value_index, names):
    quality_index = value_index + 1
    return quality_index < len(names) and (names[quality_index].startswith("Quality") or names[quality_index].startswith("Flag"))

def get_int(value):
    num = to_float(value)
    if numpy.isnan(num):
        return None
    else:
        return int(num)

def extract_data(source, index, min_vals, max_vals, replace, has_quality=False):
    if index >= 0 and is_valid_field(min_vals[index], max_vals[index]):
        return reinsert_nan(
            (numpy.nan if has_quality and get_int(d[index+1]) in [3, 4] else to_float(d[index]) for d in source),
            replace,
            length=len(source))
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

def calculate_density(
        length,
        temperature_C,
        salinity_SP,
        pressure_dbar,
        longitude,
        latitude,
        filename):
    assumed = False
    if all(x is not None for x in [temperature_C, salinity_SP, pressure_dbar]):
        salinity_SA = gsw.SA_from_SP(salinity_SP, pressure_dbar, longitude, latitude)
        density = gsw.rho(
            salinity_SA,
            gsw.CT_from_t(salinity_SA, temperature_C, pressure_dbar),
            pressure_dbar)
    else:
        logging.warning(f"Not enough data in {filename} to accurately compute density. Calculating density as though all values are 0")
        assumed = True
        density = numpy.full(length, gsw.rho([0], [0], 0)[0])
    return density, assumed

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
    density, assumed_density = calculate_density(
        len(oxygen_umol_kg),
        temperature_C,
        salinity_SP,
        pressure_dbar,
        longitude,
        latitude,
        filename)
    return numpy.fromiter(
        (o * d / oxygen_umol_per_ml * metre_cube_per_litre for o, d in zip(oxygen_umol_kg, density)),
        float,
        count=len(oxygen_umol_kg)
    ), assumed_density

def convert_salinity(
        salinity,
        units,
        filename):
    """ Converts salinity into PSU

    Arguments
    salinity - raw salinity data
    units - unit of measure for the salinity data
    filename - the name of the file the data came from (only used for logs if something goes wrong)

    Returns (oxygen in PSU, whether a computation was needed)
    """
    if salinity is None:
        return None, False
    elif units.lower() in ["psu", "pss-78"]:
        return salinity, False
    elif units.lower() in ["ppt"] or units.lower().startswith("'ppt"):
        return gsw.SP_from_SK(salinity), True
    elif units.lower() in ["umol/kg"]:
        g_per_umol = 58.44 / 1000 / 1000
        return gsw.SP_from_SR(salinity * g_per_umol), True
    else:
        warn_wrong_units("PSU", units, filename)
        return None, False

def convert_oxygen(
        oxygen,
        units,
        longitude,
        latitude,
        temperature_C,
        salinity_SP,
        pressure_dbar,
        filename):
    """ Converts oxygen concentration into mL/L

    Arguments
    oxygen - raw oxygen data
    units - unit of measure for the oxygen data
    longitude - the longitude where the oxygen data was measured
    latitude - the latitude where the oxygen data was measured
    temperature_C - the temperature (in Celcius) associated with the oxygen data
    salinity_SP - the salinity (in PSU) associated with the oxygen data
    pressure_dbar - the pressure (in dbar) associated with the oxygen data
    filename - the name of the file the data came from (only used for logs if something goes wrong)

    Returns (oxygen in mL/L, whether a computation was needed, whether density was assumed)
    """
    if oxygen is None:
        return None, False, False
    elif units.lower() in ["ml/l"]:
        return oxygen, False, False
    #                                  V this is V this when parsed with ObsFile.py
    elif units.lower() in ["umol/kg", "mmol/m", "mmol/**3"]:
        data, assumed_density = convert_umol_kg_to_mL_L(
            oxygen,
            longitude,
            latitude,
            temperature_C,
            salinity_SP,
            pressure_dbar,
            filename=filename)
        return data, True, assumed_density
    elif units.lower() in ["mg/l"]:
        oxygen_mg_per_mL = 1.429
        data = oxygen * oxygen_mg_per_mL
        return data, True, False
    else:
        warn_wrong_units("mL/L", units, filename)
        return None, False, False

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
        self.filename = os.path.basename(filename).lower()
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
            placeholder=-99.0,
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
            if datum == placeholder:
                if not once[1]:
                    logging.warning(f"Data from {filename} has value {placeholder}, which is likely a standin for NaN")
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

        if depth is None:
            if pressure is not None:
                depth = gsw.z_from_p(get_array(pressure), latitude) * -1
            else:
                logging.warning(f"{get_scalar(data.filename)} does not have depth or pressure data. Skipping")
                return

        salinity, salinity_computed = convert_salinity(
            salinity,
            None if salinity is None else salinity.units,
            filename)

        oxygen, oxygen_computed, oxygen_assumed_density = convert_oxygen(
            oxygen,
            None if oxygen is None else oxygen.units,
            longitude,
            latitude,
            temperature,
            salinity,
            pressure if pressure is not None else gsw.p_from_z(get_array(depth) * -1, latitude),
            filename)

        placeholder = -99

        if temperature is not None:
            self.temperature_data.extend(self.produce_data(
                get_array(time),
                get_array(depth),
                get_array(temperature),
                longitude,
                latitude,
                filename,
                placeholder=placeholder))

        if salinity is not None:
            self.salinity_data.extend(self.produce_data(
                get_array(time),
                get_array(depth),
                get_array(salinity),
                longitude,
                latitude,
                filename,
                placeholder=placeholder,
                computed=salinity_computed))

        if oxygen is not None:
            self.oxygen_data.extend(self.produce_data(
                get_array(time),
                get_array(depth),
                get_array(oxygen),
                longitude,
                latitude,
                filename,
                placeholder=placeholder,
                computed=oxygen_computed,
                assumed_density=oxygen_assumed_density))

    def add_data_from_shell(self, data):
        names, units, min_vals, max_vals = data.channels["Name"], data.channels["Units"], data.channels["Minimum"], data.channels["Maximum"]

        longitude, latitude = data.location["LONGITUDE"], data.location["LATITUDE"]

        time_idx = find_first(names, "Date", "DATE")
        if time_idx < 0:
            # time not included in data, just use start date
            time = numpy.full(len(data.data), data.start_dateobj)
        else:
            # time included in data
            time = extract_data(data.data, time_idx, min_vals, max_vals, data.start_dateobj)

        depth_idx = find_first(names, "Depth", "DEPTH")
        depth_pad = get_pad_value(data.channel_details, depth_idx)
        depth_data = extract_data(
            data.data,
            depth_idx,
            min_vals,
            max_vals,
            depth_pad,
            has_quality(depth_idx, names))

        temperature_idx = find_first(names, "Temperature", "TEMPERATURE")
        temperature_pad = get_pad_value(data.channel_details, temperature_idx)
        temperature_data = extract_data(
            data.data,
            temperature_idx,
            min_vals,
            max_vals,
            temperature_pad,
            has_quality(temperature_idx, names))

        salinity_idx = find_first(names, "Salinity", "SALINITY", "'Salinity", "'SALINITY")
        salinity_pad = get_pad_value(data.channel_details, salinity_idx)
        salinity_data = extract_data(
            data.data,
            salinity_idx,
            min_vals,
            max_vals,
            salinity_pad,
            has_quality(salinity_idx, names))

        oxygen_idx = find_first(names, "Oxygen", "OXYGEN")
        oxygen_pad = get_pad_value(data.channel_details, oxygen_idx)
        oxygen_data = extract_data(
            data.data,
            oxygen_idx,
            min_vals,
            max_vals,
            oxygen_pad,
            has_quality(oxygen_idx, names))

        pressure_idx = find_first(names, "Pressure", "PRESSURE")
        pressure_pad = get_pad_value(data.channel_details, pressure_idx)
        pressure_data = extract_data(
            data.data,
            pressure_idx,
            min_vals,
            max_vals,
            pressure_pad,
            has_quality(pressure_idx, names))

        if depth_data is None:
            if pressure_data is not None:
                depth_data = gsw.z_from_p(pressure_data, latitude) * -1
            else:
                logging.warning(f"{data.filename} does not have depth or pressure data. Skipping")
                return

        salinity_data, salinity_computed = convert_salinity(
            salinity_data,
            units[salinity_idx].strip(),
            data.filename)

        oxygen_data, oxygen_computed, oxygen_assumed_density = convert_oxygen(
            oxygen_data,
            units[oxygen_idx].strip(),
            longitude,
            latitude,
            temperature_data,
            salinity_data,
            pressure_data if pressure_data is not None else gsw.p_from_z(depth_data * -1, latitude),
            data.filename)

        if temperature_data is not None:
            self.temperature_data.extend(self.produce_data(
                time,
                depth_data,
                temperature_data,
                longitude,
                latitude,
                data.filename,
                placeholder=temperature_pad))

        if salinity_data is not None:
            self.salinity_data.extend(self.produce_data(
                time,
                depth_data,
                salinity_data,
                longitude,
                latitude,
                data.filename,
                placeholder=salinity_pad,
                computed=salinity_computed))

        if oxygen_data is not None:
            self.oxygen_data.extend(self.produce_data(
                time,
                depth_data,
                oxygen_data,
                longitude,
                latitude,
                data.filename,
                placeholder=oxygen_pad,
                computed=oxygen_computed,
            assumed_density=oxygen_assumed_density))
