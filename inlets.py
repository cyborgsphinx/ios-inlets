from dataclasses import dataclass
import datetime
import fnmatch
import gsw
import itertools
import json
import logging
import math
import numpy
import os
import pandas
import pickle
import re
from shapely.geometry import Point, Polygon
from typing import Dict, List
import xarray
import ios_shell.shell as ios

PICKLE_NAME = "inlets.pickle"
SHALLOW = "shallow"
MIDDLE = "middle"
DEEP = "deep"
IGNORE = "ignore"
EXCEPTIONALLY_BIG = 9.9e36


def get_length(arr):
    return arr.size if hasattr(arr, "size") else len(arr)


def is_in_bounds(val, lower, upper):
    if upper is not None:
        return lower <= val <= upper
    else:
        return lower <= val


def get_datetime(d):
    return pandas.Timestamp(d).to_pydatetime()


def reinsert_nan(data, placeholder, length=None):
    if length is None:
        length = get_length(data)
    return numpy.fromiter(
        (numpy.nan if x == placeholder or x > EXCEPTIONALLY_BIG else x for x in data),
        float,
        count=length,
    )


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


def find_column(source, name: str, *units: str) -> int:
    names = [name, "'" + name, name.upper(), "'" + name.upper()]
    potentials = [
        line for line in source if any(line.name.startswith(name) for name in names)
    ]
    units_lower = [unit.lower() for unit in units]
    with_units = [line for line in potentials if line.units.lower() in units_lower]
    if len(with_units) > 0:
        # pick first line with matching units
        return with_units[0].no - 1
    elif len(potentials) > 0:
        # pick first line with matching name
        return potentials[0].no - 1
    else:
        return -1


def to_float(source):
    if isinstance(source, float) or isinstance(source, int):
        return source
    elif isinstance(source, bytes):
        return (
            numpy.nan
            if source.strip() in [b"' '", b"n/a", b""]
            else float(source.strip().decode("utf-8"))
        )
    elif isinstance(source, str):
        return (
            numpy.nan if source.strip() in ["' '", "n/a", ""] else float(source.strip())
        )
    else:
        raise ValueError(f"to_float called on {source}")


def find_all(source, attrs):
    out = []
    for attr in attrs:
        if hasattr(source, attr):
            out.append(getattr(source, attr))
    return out


def find_data(source, names, units):
    potentials = find_all(source, names)
    with_units = [
        column
        for column in potentials
        if isinstance(column, float) or column.units.lower() in units
    ]
    if len(with_units) > 0:
        return with_units[0]
    elif len(potentials) > 0:
        return potentials[0]
    else:
        return None


def find_temperature_data(data):
    temperature_names = [
        "TEMPRTN1",
        "TEMPST01",
        "TEMPPR01",
        "TEMPPR03",
        "TEMPS901",
        "TEMPS601",
    ]
    temperature_units = ["C", "deg C", "degrees C"]
    return find_data(data, temperature_names, temperature_units)


def find_salinity_data(data):
    salinity_names = [
        "PSLTZZ01",
        "ODSDM021",
        "SSALST01",
        "PSALST01",
        "PSALBST1",
        "sea_water_practical_salinity",
    ]
    salinity_units = ["PSU", "PSS-78"]
    return find_data(data, salinity_names, salinity_units)


def find_oxygen_data(data):
    oxygen_names = ["DOXYZZ01", "DOXMZZ01"]
    oxygen_units = ["mL/L"]
    return find_data(data, oxygen_names, oxygen_units)


def find_depth_data(data):
    depth_names = ["depth", "instrument_depth", "PPSAADCP"]
    depth_units = ["m", "metres"]
    return find_data(data, depth_names, depth_units)


def find_pressure_data(data):
    pressure_names = ["PRESPR01", "sea_water_pressure"]
    pressure_units = ["dbar", "decibar", "decibars"]
    return find_data(data, pressure_names, pressure_units)


def extend_arr(arr, length):
    arr_length = get_length(arr)
    return numpy.full(length, arr.item()) if arr_length == 1 else arr


def get_pad_value(info, index):
    if index < 0 or info is None or len(info) == 0:
        return None
    return to_float(info[index].pad)


def has_quality(value_index, names):
    quality_index = value_index + 1
    return quality_index < len(names) and (
        names[quality_index].startswith("Quality")
        or names[quality_index].startswith("Flag")
    )


def is_acceptable_quality(quality_value):
    # 2 is "inconsistent with climatology" in the vast majority of observed cases
    bad_qualities = ["2", "3", "4"]
    return not any(value in quality_value for value in bad_qualities)


def get_int(value):
    num = to_float(value)
    if numpy.isnan(num):
        return None
    else:
        return int(num)


def extract_data(source, index, replace, has_quality=False):
    if index < 0:
        return None
    return reinsert_nan(
        (
            numpy.nan
            if has_quality and not is_acceptable_quality(str(row[index + 1]))
            else to_float(row[index])
            for row in source
        ),
        replace,
        length=len(source),
    )


def warn_unknown_variable(data, var):
    # check if there is a potential variable based on the broader name
    var_list = []
    for key in data.keys():
        if re.search(var, getattr(data[key], "long_name", "").lower()):
            var_list.append(key)
    if len(var_list) != 0:
        logging.warning(
            f"{get_scalar(data.filename)} has unknown {var} variable. Possible values: {var_list}"
        )


def warn_wrong_units(expected, actual, filename):
    logging.warning(
        f"Cowardly refusing to perform the conversion from {actual} to {expected} in {filename}"
    )


def calculate_density(
    length, temperature_C, salinity_SP, pressure_dbar, longitude, latitude, filename
):
    assumed = False
    if all(x is not None for x in [temperature_C, salinity_SP, pressure_dbar]):
        salinity_SA = gsw.SA_from_SP(salinity_SP, pressure_dbar, longitude, latitude)
        density = gsw.rho(
            salinity_SA,
            gsw.CT_from_t(salinity_SA, temperature_C, pressure_dbar),
            pressure_dbar,
        )
    else:
        logging.warning(
            f"Not enough data in {filename} to accurately compute density. Calculating density as though all values are 0"
        )
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
    filename="unknown file",
):
    oxygen_umol_per_ml = 44.661
    metre_cube_per_litre = 0.001
    density, assumed_density = calculate_density(
        len(oxygen_umol_kg),
        temperature_C,
        salinity_SP,
        pressure_dbar,
        longitude,
        latitude,
        filename,
    )
    return (
        numpy.fromiter(
            (
                o * d / oxygen_umol_per_ml * metre_cube_per_litre
                for o, d in zip(oxygen_umol_kg, density)
            ),
            float,
            count=len(oxygen_umol_kg),
        ),
        assumed_density,
    )


def convert_percent_to_mL_L(
    oxygen_percent,
    temperature_C,
    salinity_SP,
    filename="no filename given",
):
    if temperature_C is not None and salinity_SP is not None:
        # function from en.wikipedia.org/wiki/Oxygenation_(environmental)
        kelven_offset = 273.15
        temperature_K = [t + kelven_offset for t in temperature_C]
        A1 = -173.4292
        A2 = 249.6339
        A3 = 143.3483
        A4 = -21.8492
        B1 = -0.033096
        B2 = 0.014259
        B3 = -0.001700
        return numpy.fromiter(
            (
                (o / 100)
                * math.exp(
                    A1
                    + (A2 * 100 / t)
                    + (A3 * math.log(t / 100))
                    + (A4 * t / 100)
                    + (s * (B1 + (B2 * t / 100) + (B3 * ((t / 100) ** 2))))
                )
                for o, t, s in zip(oxygen_percent, temperature_K, salinity_SP)
            ),
            float,
            count=len(oxygen_percent),
        )
    else:
        logging.warning(
            f"Not enough data from {filename} to convert oxygen from % to mL/L. Ignoring file"
        )
        return None


def convert_salinity(salinity, units, filename):
    """Converts salinity into PSU

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
    filename,
):
    """Converts oxygen concentration into mL/L

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
    elif units.lower() in ["umol/kg", "mmol/m", "mmol/m**3"]:
        data, assumed_density = convert_umol_kg_to_mL_L(
            oxygen,
            longitude,
            latitude,
            temperature_C,
            salinity_SP,
            pressure_dbar,
            filename=filename,
        )
        return data, True, assumed_density
    elif units.lower() in ["mg/l"]:
        oxygen_mg_per_mL = 1.429
        data = oxygen * oxygen_mg_per_mL
        return data, True, False
    elif units in ["%"]:
        data = convert_percent_to_mL_L(
            oxygen, temperature_C, salinity_SP, filename=filename
        )
        return data, False, False
    else:
        warn_wrong_units("mL/L", units, filename)
        return None, False, False


def get_data(col, bucket, before=None):
    data = [
        [datum.time, datum.datum]
        for datum in col
        if datum.bucket == bucket and not numpy.isnan(datum.datum)
    ]
    if before is not None:
        data = [[t, d] for t, d in data if t.year < before.year]
    return zip(*data) if len(data) > 0 else [[], []]


@dataclass
class InletData(object):
    time: datetime.datetime
    bucket: str
    depth: float
    datum: float
    longitude: float
    latitude: float
    filename: str
    computed: bool = False
    assumed_density: bool = False


def get_outliers(col, bucket="all", before=None):
    data = [
        datum
        for datum in col
        if (bucket == "all" or datum.bucket == bucket) and not numpy.isnan(datum.datum)
    ]
    if before is not None:
        data = [datum for datum in data if datum.time.year < before.year]
    return [datum for datum in data if datum.datum < 3]


class Inlet(object):
    def __init__(
        self,
        name: str,
        polygon: Polygon,
        boundaries: List[int],
        limits: Dict[str, List[float]],
    ):
        self.name = name
        self.shallow_bounds = (boundaries[0], boundaries[1])
        self.middle_bounds = (boundaries[1], boundaries[2])
        self.deep_bounds = (boundaries[2], None)
        self.temperature_data = []
        self.salinity_data = []
        self.oxygen_data = []
        self.polygon = polygon
        self.limits = limits
        self.used_files = set()

    def get_temperature_data(self, bucket, before=None):
        return get_data(self.temperature_data, bucket, before)

    def get_salinity_data(self, bucket, before=None):
        return get_data(self.salinity_data, bucket, before)

    def get_oxygen_data(self, bucket, before=None):
        return get_data(self.oxygen_data, bucket, before)

    def get_temperature_outliers(self, bucket="all", before=None):
        return get_outliers(self.temperature_data, bucket, before)

    def get_salinity_outliers(self, bucket="all", before=None):
        return get_outliers(self.salinity_data, bucket, before)

    def get_oxygen_outliers(self, bucket="all", before=None):
        return get_outliers(self.oxygen_data, bucket, before)

    def has_temperature_data(self):
        return len(self.temperature_data) > 0

    def has_salinity_data(self):
        return len(self.salinity_data) > 0

    def has_oxygen_data(self):
        return len(self.oxygen_data) > 0

    def has_data_from(self, file_name):
        return os.path.basename(file_name).lower() in self.used_files

    def get_station_data(self, before=None):
        temperature_data = (
            filter(lambda x: x.time.year < before.year, self.temperature_data)
            if before is not None
            else self.temperature_data
        )
        salinity_data = (
            filter(lambda x: x.time.year < before.year, self.salinity_data)
            if before is not None
            else self.salinity_data
        )
        oxygen_data = (
            filter(lambda x: x.time.year < before.year, self.oxygen_data)
            if before is not None
            else self.oxygen_data
        )
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
        assumed_density=False,
    ):
        length = get_length(data)
        times = extend_arr(times, length)
        depths = extend_arr(depths, length)
        if (get_length(times) != get_length(data)) or (
            get_length(depths) != get_length(data)
        ):
            logging.warning(
                f"Data from {filename} contains times, depths, and data of different lengths"
            )

        out = []
        once = [False] * 3
        for t, d, datum in zip(times, depths, data):
            # Some data, particularly salinity data, seems to be the result of performing calculations on NaN values.
            # This data is consistently showing up as 9.96921e+36, which may relate to the "Fill Value" in creating netCDF files.
            # In any case, it appears to be as invalid as NaN, so it's being filtered out accordingly
            if datum > EXCEPTIONALLY_BIG or d > EXCEPTIONALLY_BIG:
                if not once[0]:
                    logging.warning(
                        f"Data from {filename} is larger than 9.9e+36, it may have been calulated poorly"
                    )
                    once[0] = True
                continue
            if datum == placeholder:
                if not once[1]:
                    logging.warning(
                        f"Data from {filename} has value {placeholder}, which is likely a standin for NaN"
                    )
                    once[1] = True
                datum = numpy.nan
            if self.is_shallow(d):
                category = SHALLOW
            elif self.is_middle(d):
                category = MIDDLE
            elif self.is_deep(d):
                category = DEEP
            else:
                category = IGNORE
            out.append(
                InletData(
                    get_datetime(t),
                    category,
                    d,
                    datum,
                    longitude,
                    latitude,
                    filename,
                    computed=computed,
                    assumed_density=assumed_density,
                )
            )
        if len(out) == 0:
            logging.warning(f"Data from {filename} not used")
        else:
            self.used_files.add(os.path.basename(filename).lower())
        return out

    def add_data_from_netcdf(self, data):
        time, longitude, latitude, filename = (
            get_array(data.time),
            get_scalar(data.longitude),
            get_scalar(data.latitude),
            get_scalar(data.filename),
        )

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
                logging.warning(
                    f"{get_scalar(data.filename)} does not have depth or pressure data. Treating depth as NaN"
                )
                depth = numpy.nan

        salinity, salinity_computed = convert_salinity(
            salinity,
            None if salinity is None or isinstance(salinity, float) else salinity.units,
            filename,
        )

        oxygen, oxygen_computed, oxygen_assumed_density = convert_oxygen(
            oxygen,
            None if oxygen is None or isinstance(oxygen, float) else oxygen.units,
            longitude,
            latitude,
            temperature,
            salinity,
            pressure
            if pressure is not None
            else gsw.p_from_z(get_array(depth) * -1, latitude),
            filename,
        )

        placeholder = -99

        if temperature is not None:
            self.temperature_data.extend(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(temperature),
                    longitude,
                    latitude,
                    filename,
                    placeholder=placeholder,
                )
            )

        if salinity is not None:
            self.salinity_data.extend(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(salinity),
                    longitude,
                    latitude,
                    filename,
                    placeholder=placeholder,
                    computed=salinity_computed,
                )
            )

        if oxygen is not None:
            self.oxygen_data.extend(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(oxygen),
                    longitude,
                    latitude,
                    filename,
                    placeholder=placeholder,
                    computed=oxygen_computed,
                    assumed_density=oxygen_assumed_density,
                )
            )

    def add_data_from_shell(self, data):
        channels = data.file.channels
        channel_details = data.file.channel_details
        names = [channel.name for channel in channels]
        units = [channel.units for channel in channels]

        longitude, latitude = data.location.longitude, data.location.latitude

        date_idx = find_column(channels, "Date")
        if date_idx < 0:
            # time not included in data, just use start date
            time = numpy.full(len(data.data), data.get_time())
        else:
            time_idx = find_column(channels, "Time")
            if time_idx < 0:
                # only date included in data
                time = [d[date_idx] for d in data.data]
            else:
                dates = [d[date_idx] for d in data.data]
                times = [d[time_idx] for d in data.data]
                time = [
                    datetime.datetime.combine(d, t, tzinfo=data.get_time().tzinfo)
                    for d, t in zip(dates, times)
                ]

        depth_idx = find_column(channels, "Depth", "m", "metre")
        depth_pad = get_pad_value(channel_details, depth_idx)
        if depth_pad is None or numpy.isnan(depth_pad):
            depth_pad = -99
        depth_data = extract_data(
            data.data, depth_idx, depth_pad, has_quality(depth_idx, names)
        )

        temperature_idx = find_column(channels, "Temperature", "C", "'deg C'")
        temperature_pad = get_pad_value(channel_details, temperature_idx)
        if temperature_pad is None or numpy.isnan(temperature_pad):
            temperature_pad = -99
        temperature_data = extract_data(
            data.data,
            temperature_idx,
            temperature_pad,
            has_quality(temperature_idx, names),
        )

        salinity_idx = find_column(channels, "Salinity", "PSU", "PSS-78")
        salinity_pad = get_pad_value(channel_details, salinity_idx)
        if salinity_pad is None or numpy.isnan(salinity_pad):
            salinity_pad = -99
        salinity_data = extract_data(
            data.data, salinity_idx, salinity_pad, has_quality(salinity_idx, names)
        )

        oxygen_idx = find_column(channels, "Oxygen", "mL/L")
        oxygen_pad = get_pad_value(channel_details, oxygen_idx)
        if oxygen_pad is None or numpy.isnan(oxygen_pad):
            oxygen_pad = -99
        oxygen_data = extract_data(
            data.data, oxygen_idx, oxygen_pad, has_quality(oxygen_idx, names)
        )

        pressure_idx = find_column(channels, "Pressure", "dbar", "decibar")
        pressure_pad = get_pad_value(channel_details, pressure_idx)
        if pressure_pad is None or numpy.isnan(pressure_pad):
            pressure_pad = -99
        pressure_data = extract_data(
            data.data, pressure_idx, pressure_pad, has_quality(pressure_idx, names)
        )

        if (
            depth_data is None
            and data.instrument is not None
            and not numpy.isnan(data.instrument.depth)
        ):
            depth_data = numpy.full(1, float(data.instrument.raw["depth"]))
        elif depth_data is None:
            if pressure_data is not None:
                depth_data = gsw.z_from_p(pressure_data, latitude) * -1
            else:
                logging.warning(
                    f"{data.filename} does not have depth or pressure data. Skipping"
                )
                return
        elif pressure_data is None:
            # depth_data is not None in this case
            pressure_data = gsw.p_from_z(depth_data * -1, latitude)

        salinity_data, salinity_computed = convert_salinity(
            salinity_data, units[salinity_idx].strip(), data.filename
        )

        oxygen_data, oxygen_computed, oxygen_assumed_density = convert_oxygen(
            oxygen_data,
            units[oxygen_idx].strip(),
            longitude,
            latitude,
            temperature_data,
            salinity_data,
            pressure_data
            if pressure_data is not None
            else gsw.p_from_z(depth_data * -1, latitude),
            data.filename,
        )

        if temperature_data is not None:
            self.temperature_data.extend(
                self.produce_data(
                    time,
                    depth_data,
                    temperature_data,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=temperature_pad,
                )
            )

        if salinity_data is not None:
            self.salinity_data.extend(
                self.produce_data(
                    time,
                    depth_data,
                    salinity_data,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=salinity_pad,
                    computed=salinity_computed,
                )
            )

        if oxygen_data is not None:
            self.oxygen_data.extend(
                self.produce_data(
                    time,
                    depth_data,
                    oxygen_data,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=oxygen_pad,
                    computed=oxygen_computed,
                    assumed_density=oxygen_assumed_density,
                )
            )


def get_inlets(data_dir, from_saved=False, skip_netcdf=False):
    inlet_list = []
    if from_saved:
        with open(PICKLE_NAME, mode="rb") as f:
            inlet_list = pickle.load(f)
    else:
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
                inlet_list.append(Inlet(name, polygon, boundaries, limits))

        if not skip_netcdf:
            for root, _, files in os.walk(os.path.join(data_dir, "netCDF_Data")):
                for item in fnmatch.filter(files, "*.nc"):
                    file_name = os.path.join(root, item)
                    data = xarray.open_dataset(file_name)
                    for inlet in inlet_list:
                        if inlet.contains(data):
                            try:
                                inlet.add_data_from_netcdf(data)
                            except:
                                logging.exception(f"Exception occurred in {file_name}")
                                raise

        shell_exts = ["bot", "che", "ctd", "ubc", "med", "xbt", "adcp", "cur"]
        # make a list of all elements in shell_exts followed by their str.upper() versions
        exts = [
            item
            for sublist in [[ext, ext.upper()] for ext in shell_exts]
            for item in sublist
        ]
        for root, dirs, files in os.walk(data_dir):
            for ext in exts:
                for item in fnmatch.filter(files, "*." + ext):
                    file_name = os.path.join(root, item)
                    try:
                        shell = ios.ShellFile.fromfile(file_name, process_data=False)
                    except Exception as e:
                        logging.exception(f"Failed to read {file_name}: {e}")
                        continue
                    for inlet in inlet_list:
                        if inlet.contains(shell.get_location()):
                            # use item instead of file_name because the netcdf files don't store path information
                            # they also do not store the .nc extension, so this should be reasonable
                            if not inlet.has_data_from(item.lower()):
                                try:
                                    shell.process_data()
                                    inlet.add_data_from_shell(shell)
                                except Exception as e:
                                    logging.exception(
                                        f"Exception occurred in {file_name}: {e}"
                                    )
                                    continue
            if "HISTORY" in dirs:
                dirs.remove("HISTORY")
        with open(PICKLE_NAME, mode="wb") as f:
            pickle.dump(inlet_list, f)
    return inlet_list
