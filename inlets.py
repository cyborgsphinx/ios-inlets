import convert
import csv
import datetime
import fnmatch
import gsw
import inlet_data
import itertools
import json
import logging
import math
import numpy
import os
import pandas
import re
from shapely.geometry import Point, Polygon
from typing import Dict, List
import xarray
import ios_shell.shell as ios
import erddap
from enum import Enum

EXCEPTIONALLY_BIG = 9.9e36

class Category(Enum):
    IGNORE = 1
    SURFACE = 2
    SHALLOW = 3
    USED_SURFACE = 4
    DEEP = 5
    DEEPER = 6
    DEEPEST = 7
    USED_DEEP = 8
    ALL = 9


def get_length(arr):
    if arr is None:
        return 0
    return arr.size if hasattr(arr, "size") else len(arr)


def is_in_bounds(val, lower, upper):
    if upper is not None:
        return lower <= val <= upper
    else:
        return lower <= val


def get_datetime(d):
    return pandas.to_datetime(d).to_pydatetime()


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
    name_lower = name.lower()
    potentials = [line for line in source if name_lower in line.name.lower()]
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
    depth_names = ["depth", "depth_nominal", "instrument_depth", "PPSAADCP"]
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
    bad_qualities = [2, 3, 4]
    return quality_value not in bad_qualities


def extract_data(source, index, replace):
    if index < 0:
        return None
    return reinsert_nan(
        (to_float(row[index]) for row in source),
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


def get_data(col, before=None, do_average=False):
    data = [
        [datum.time, datum.value]
        for datum in col
        if is_acceptable_quality(datum.quality)
    ]
    if before is not None:
        data = [[t, d] for t, d in data if t.year < before.year]
    if do_average:
        data_dict = {}
        for t, d in data:
            date = datetime.date(t.year, t.month, 1)
            if date not in data_dict:
                data_dict[date] = {"total": 0, "count": 0}
            data_dict[date]["total"] += d
            data_dict[date]["count"] += 1
        data = [[key, elem["total"] / elem["count"]] for key, elem in data_dict.items()]
    return zip(*data) if len(data) > 0 else [[], []]


def hakai_quality(quality):
    del quality
    # assume all qualities are good for now
    return 1


class Inlet(object):
    def __init__(
        self,
        name: str,
        area: str,
        polygon: Polygon,
        boundaries: List[int],
        limits: Dict[str, List[float]],
        clear_old_data: bool = False,
        db_name=None,
        shallow: List[int] = [0, 30, 100],
    ):
        self.name = name
        self.area = area
        self.deep_bounds = (boundaries[0], boundaries[1])
        self.deeper_bounds = (boundaries[1], boundaries[2])
        self.deepest_bounds = (
            boundaries[2],
            boundaries[3] if len(boundaries) > 3 else None,
        )
        self.polygon = polygon
        self.limits = limits
        self.used_files = set()
        if db_name is not None:
            self.data = inlet_data.InletDb(name, clear_old_data, db_name)
        else:
            self.data = inlet_data.InletDb(name, clear_old_data)
        self.surface_bounds = (shallow[0], shallow[1])
        if len(shallow) > 2:
            self.shallow_bounds = (shallow[1], shallow[2])
        else:
            self.shallow_bounds = None

    def __bucket_to_bounds(self, bucket: Category):
        return (
            (None, None)
            if bucket == Category.ALL
            else self.surface_bounds
            if bucket == Category.SURFACE
            or bucket == Category.USED_SURFACE
            and self.shallow_bounds is None
            else self.shallow_bounds
            if bucket == Category.SHALLOW
            else self.deep_bounds
            if bucket == Category.DEEP
            else self.deeper_bounds
            if bucket == Category.DEEPER
            else self.deepest_bounds
            if bucket == Category.DEEPEST
            else (self.shallow_bounds[1], self.deep_bounds[0])
            if bucket == Category.IGNORE
            else (self.surface_bounds[0], self.shallow_bounds[1])
            if bucket == Category.USED_SURFACE
            else (self.deep_bounds[0], self.deepest_bounds[1])
            if bucket == Category.USED_DEEP
            else None
        )

    def get_temperature_data(self, bucket: str, before=None, do_average=False):
        if bucket == Category.SHALLOW and self.shallow_bounds is None:
            return [[], []]
        bounds = self.__bucket_to_bounds(bucket)
        return get_data(
            self.data.get_temperature_data(bounds, average=do_average),
            before,
            do_average,
        )

    def get_salinity_data(self, bucket: str, before=None, do_average=False):
        if bucket == Category.SHALLOW and self.shallow_bounds is None:
            return [[], []]
        bounds = self.__bucket_to_bounds(bucket)
        return get_data(
            self.data.get_salinity_data(bounds, average=do_average), before, do_average
        )

    def get_oxygen_data(self, bucket: str, before=None, do_average=False):
        if bucket == Category.SHALLOW and self.shallow_bounds is None:
            return [[], []]
        bounds = self.__bucket_to_bounds(bucket)
        return get_data(
            self.data.get_oxygen_data(bounds, average=do_average), before, do_average
        )

    def has_temperature_data(self):
        bounds = self.__bucket_to_bounds(Category.ALL)
        return len(self.data.get_temperature_data(bounds)) > 0

    def has_salinity_data(self):
        bounds = self.__bucket_to_bounds(Category.ALL)
        return len(self.data.get_salinity_data(bounds)) > 0

    def has_oxygen_data(self):
        bounds = self.__bucket_to_bounds(Category.ALL)
        return len(self.data.get_oxygen_data(bounds)) > 0

    def has_data_from(self, file_name):
        return os.path.basename(file_name).lower() in self.used_files

    def get_station_data(self, before=None):
        bounds = self.__bucket_to_bounds(Category.ALL)
        data = self.data.get_temperature_data(bounds)
        temperature_data = (
            filter(lambda x: x.time.year < before.year, data)
            if before is not None
            else data
        )
        data = self.data.get_salinity_data(bounds)
        salinity_data = (
            filter(lambda x: x.time.year < before.year, data)
            if before is not None
            else data
        )
        data = self.data.get_oxygen_data(bounds)
        oxygen_data = (
            filter(lambda x: x.time.year < before.year, data)
            if before is not None
            else data
        )
        stations = {}
        for datum in itertools.chain(temperature_data, salinity_data, oxygen_data):
            year = datum.time.year
            if year not in stations:
                stations[year] = set()
            stations[year].add(datum.source)
        return stations

    def contains(self, latitude=None, longitude=None):
        if longitude is None:
            logging.warning("data does not contain longitude information")
            return False

        if latitude is None:
            logging.warning("data does not contain latitude information")
            return False

        return self.polygon.contains(Point(longitude, latitude))

    def bounding_box(self):
        min_lon, min_lat, max_lon, max_lat = self.polygon.bounds
        return {
            "min_lon": min_lon,
            "max_lon": max_lon,
            "min_lat": min_lat,
            "max_lat": max_lat,
        }

    def is_surface(self, depth):
        return is_in_bounds(depth, *self.surface_bounds)

    def is_shallow(self, depth):
        return self.shallow_bounds is not None and is_in_bounds(
            depth, *self.shallow_bounds
        )

    def is_deep(self, depth):
        return is_in_bounds(depth, *self.deep_bounds)

    def is_deeper(self, depth):
        return is_in_bounds(depth, *self.deeper_bounds)

    def is_deepest(self, depth):
        return is_in_bounds(depth, *self.deepest_bounds)

    def produce_data(
        self,
        times,
        depths,
        data,
        quality,
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
        once = [False] * 2
        warn_unused = True
        for t, d, datum, q in zip(times, depths, data, quality):
            if math.isnan(datum):
                # no warning since NaN data is incredibly common
                # if a file winds up with no data because all the data was NaN, don't warn that it wasn't used
                warn_unused = False
                continue
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
            if datum == placeholder or int(datum) == int(placeholder):
                if not once[1]:
                    logging.warning(
                        f"Data from {filename} has value {datum}, which is likely a standin for NaN"
                    )
                    once[1] = True
                continue
            t = get_datetime(t)
            if t.replace(tzinfo=None) > datetime.datetime.now():
                logging.warning(f"Data from {filename} is from the future: {t}")
                continue
            out.append(
                inlet_data.InletData(
                    t,
                    d,
                    datum,
                    q if q is not None and math.isfinite(q) else 0,
                    longitude,
                    latitude,
                    filename,
                    computed=computed,
                    assumed_density=assumed_density,
                )
            )
        if len(out) == 0:
            if warn_unused:
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

        salinity, salinity_computed = convert.convert_salinity(
            salinity,
            None if salinity is None or isinstance(salinity, float) else salinity.units,
            filename,
        )

        oxygen, oxygen_computed, oxygen_assumed_density = convert.convert_oxygen(
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
        assumed_quality = 1  # assume good for netCDF data

        if temperature is not None:
            self.data.add_temperature_data(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(temperature),
                    numpy.full(get_length(temperature), assumed_quality),
                    longitude,
                    latitude,
                    filename,
                    placeholder=placeholder,
                )
            )

        if salinity is not None:
            self.data.add_salinity_data(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(salinity),
                    numpy.full(get_length(salinity), assumed_quality),
                    longitude,
                    latitude,
                    filename,
                    placeholder=placeholder,
                    computed=salinity_computed,
                )
            )

        if oxygen is not None:
            self.data.add_oxygen_data(
                self.produce_data(
                    get_array(time),
                    get_array(depth),
                    get_array(oxygen),
                    numpy.full(get_length(oxygen), assumed_quality),
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
        depth_data = extract_data(data.data, depth_idx, depth_pad)

        temperature_idx = find_column(channels, "Temperature", "C", "'deg C'")
        temperature_pad = get_pad_value(channel_details, temperature_idx)
        if temperature_pad is None or numpy.isnan(temperature_pad):
            temperature_pad = -99
        temperature_data = extract_data(data.data, temperature_idx, temperature_pad)
        temperature_quality = [0] * get_length(temperature_data)
        if has_quality(temperature_idx, names):
            temperature_quality = extract_data(data.data, temperature_idx + 1, None)

        salinity_idx = find_column(channels, "Salinity", "PSU", "PSS-78")
        salinity_pad = get_pad_value(channel_details, salinity_idx)
        if salinity_pad is None or numpy.isnan(salinity_pad):
            salinity_pad = -99
        salinity_data = extract_data(data.data, salinity_idx, salinity_pad)
        salinity_quality = [0] * get_length(salinity_data)
        if has_quality(salinity_idx, names):
            salinity_quality = extract_data(data.data, salinity_idx + 1, None)

        oxygen_idx = find_column(channels, "Oxygen", "mL/L")
        oxygen_pad = get_pad_value(channel_details, oxygen_idx)
        if oxygen_pad is None or numpy.isnan(oxygen_pad):
            oxygen_pad = -99
        oxygen_data = extract_data(data.data, oxygen_idx, oxygen_pad)
        oxygen_quality = [0] * get_length(oxygen_data)
        if has_quality(oxygen_idx, names):
            oxygen_quality = extract_data(data.data, oxygen_idx + 1, None)

        pressure_idx = find_column(channels, "Pressure", "dbar", "decibar")
        pressure_pad = get_pad_value(channel_details, pressure_idx)
        if pressure_pad is None or numpy.isnan(pressure_pad):
            pressure_pad = -99
        pressure_data = extract_data(data.data, pressure_idx, pressure_pad)

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

        salinity_data, salinity_computed = convert.convert_salinity(
            salinity_data, units[salinity_idx].strip(), data.filename
        )

        oxygen_data, oxygen_computed, oxygen_assumed_density = convert.convert_oxygen(
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
            self.data.add_temperature_data(
                self.produce_data(
                    time,
                    depth_data,
                    temperature_data,
                    temperature_quality,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=temperature_pad,
                )
            )

        if salinity_data is not None:
            self.data.add_salinity_data(
                self.produce_data(
                    time,
                    depth_data,
                    salinity_data,
                    salinity_quality,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=salinity_pad,
                    computed=salinity_computed,
                )
            )

        if oxygen_data is not None:
            self.data.add_oxygen_data(
                self.produce_data(
                    time,
                    depth_data,
                    oxygen_data,
                    oxygen_quality,
                    longitude,
                    latitude,
                    data.filename,
                    placeholder=oxygen_pad,
                    computed=oxygen_computed,
                    assumed_density=oxygen_assumed_density,
                )
            )

    def add_data_from_csv(self, data, filename):
        time = pandas.to_datetime(data["Measurement time"])
        longitude = data["Longitude"]
        latitude = data["Latitude"]
        depth = data["Depth (m)"]
        temperature = data["Temperature (deg C)"]
        temperature_flag = data["Temperature flag"].map(hakai_quality)
        oxygen_ml_l = data["Dissolved O2 (mL/L)"]
        oxygen_flag = data["Dissolved O2 (mL/L) flag"].map(hakai_quality)
        salinity = data["Salinity (PSU)"]
        salinity_flag = data["Salinity flag"].map(hakai_quality)

        temperature_index = temperature.notna()
        salinity_index = salinity.notna()
        oxygen_index = oxygen_ml_l.notna()

        self.data.add_temperature_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                )
                for t, d, v, q, lon, lat in zip(
                    time[temperature_index],
                    depth[temperature_index],
                    temperature[temperature_index],
                    temperature_flag[temperature_index],
                    longitude[temperature_index],
                    latitude[temperature_index],
                )
            ]
        )

        self.data.add_salinity_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                )
                for t, d, v, q, lon, lat in zip(
                    time[salinity_index],
                    depth[salinity_index],
                    salinity[salinity_index],
                    salinity_flag[salinity_index],
                    longitude[salinity_index],
                    latitude[salinity_index],
                )
            ]
        )

        self.data.add_oxygen_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                )
                for t, d, v, q, lon, lat in zip(
                    time[oxygen_index],
                    depth[oxygen_index],
                    oxygen_ml_l[oxygen_index],
                    oxygen_flag[oxygen_index],
                    longitude[oxygen_index],
                    latitude[oxygen_index],
                )
            ]
        )

    def add_data_from_erddap(self, data):
        depth_index = data["depth"].map(math.isfinite)
        temperature_index = (
            data["aggregated_temperature"].map(math.isfinite) & depth_index
        )
        salinity_index = data["aggregated_salinity"].map(math.isfinite) & depth_index
        oxygen_index = data["aggregated_oxygen"].map(math.isfinite) & depth_index

        self.data.add_temperature_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                    computed=computed,
                    assumed_density=assumed,
                )
                for t, d, v, q, lon, lat, filename, computed, assumed in zip(
                    data.loc[temperature_index, "time"].map(get_datetime),
                    data.loc[temperature_index, "depth"],
                    data.loc[temperature_index, "aggregated_temperature"],
                    data.loc[temperature_index, "aggregated_temperature_quality"],
                    data.loc[temperature_index, "longitude"],
                    data.loc[temperature_index, "latitude"],
                    data.loc[temperature_index, "source"],
                    data.loc[temperature_index, "aggregated_temperature_metadata"].map(
                        lambda x: x > 1
                    ),
                    data.loc[temperature_index, "aggregated_temperature_metadata"].map(
                        lambda x: x > 2
                    ),
                )
            ]
        )

        self.data.add_salinity_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                    computed=computed,
                    assumed_density=assumed,
                )
                for t, d, v, q, lon, lat, filename, computed, assumed in zip(
                    data.loc[salinity_index, "time"].map(get_datetime),
                    data.loc[salinity_index, "depth"],
                    data.loc[salinity_index, "aggregated_salinity"],
                    data.loc[salinity_index, "aggregated_salinity_quality"],
                    data.loc[salinity_index, "longitude"],
                    data.loc[salinity_index, "latitude"],
                    data.loc[salinity_index, "source"],
                    data.loc[salinity_index, "aggregated_salinity_metadata"].map(
                        lambda x: x > 1
                    ),
                    data.loc[salinity_index, "aggregated_salinity_metadata"].map(
                        lambda x: x > 2
                    ),
                )
            ]
        )

        self.data.add_oxygen_data(
            [
                inlet_data.InletData(
                    time=t,
                    depth=d,
                    value=v,
                    quality=q,
                    longitude=lon,
                    latitude=lat,
                    source=filename,
                    computed=computed,
                    assumed_density=assumed,
                )
                for t, d, v, q, lon, lat, filename, computed, assumed in zip(
                    data.loc[oxygen_index, "time"].map(get_datetime),
                    data.loc[oxygen_index, "depth"],
                    data.loc[oxygen_index, "aggregated_oxygen"],
                    data.loc[oxygen_index, "aggregated_oxygen_quality"],
                    data.loc[oxygen_index, "longitude"],
                    data.loc[oxygen_index, "latitude"],
                    data.loc[oxygen_index, "source"],
                    data.loc[oxygen_index, "aggregated_oxygen_metadata"].map(
                        lambda x: x > 1
                    ),
                    data.loc[oxygen_index, "aggregated_oxygen_metadata"].map(
                        lambda x: x > 2
                    ),
                )
            ]
        )


def get_inlets(
    data_dir,
    from_saved=False,
    from_netcdf=False,
    from_erddap=False,
    from_csv=False,
    inlet_names=[],
    drop_names=[],
    keep_names=[],
    geojson_file="inlets.geojson",
) -> List[Inlet]:
    inlet_list = []
    with open(geojson_file) as f:
        contents = json.load(f)["features"]
        for content in contents:
            name = content["properties"]["name"]
            if len(keep_names) > 0 and not all(
                name_part in name for name_part in keep_names
            ):
                continue
            if len(inlet_names) > 0 and not any(
                name_part in name for name_part in inlet_names
            ):
                continue
            if len(drop_names) > 0 and any(
                name_part in name for name_part in drop_names
            ):
                continue
            area = content["properties"]["area"]
            boundaries = content["properties"]["boundaries"]
            limits = (
                content["properties"]["limits"]
                if "limits" in content["properties"]
                else {}
            )
            polygon = Polygon(content["geometry"]["coordinates"][0])
            if "shallow boundaries" in content["properties"]:
                inlet_list.append(
                    Inlet(
                        name,
                        area,
                        polygon,
                        boundaries,
                        limits,
                        clear_old_data=not from_saved,
                        shallow=content["properties"]["shallow boundaries"],
                    )
                )
            else:
                inlet_list.append(
                    Inlet(
                        name,
                        area,
                        polygon,
                        boundaries,
                        limits,
                        clear_old_data=not from_saved,
                    )
                )
    if not from_saved:
        if from_erddap:
            for inlet in inlet_list:
                for data_frame in erddap.pull_data_for(inlet):
                    inlet.add_data_from_erddap(data_frame)

        if from_netcdf:
            for root, _, files in os.walk(os.path.join(data_dir, "netCDF_Data")):
                for item in fnmatch.filter(files, "*.nc"):
                    file_name = os.path.join(root, item)
                    data = xarray.open_dataset(file_name)
                    for inlet in inlet_list:
                        if inlet.contains(
                            latitude=data.latitude, longitude=data.longitude
                        ):
                            try:
                                inlet.add_data_from_netcdf(data)
                            except:
                                logging.exception(f"Exception occurred in {file_name}")
                                raise

        shell_exts = ["bot", "che", "ctd", "ubc", "med", "xbt", "cur"]
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
                    except Exception:
                        logging.exception(f"Error encountered reading {file_name}")
                        continue
                    for inlet in inlet_list:
                        if inlet.contains(**shell.get_location()):
                            # Use item instead of file_name because the netcdf files don't store
                            # path information. They also do not store the .nc extension, so this
                            # should be reasonable
                            if not inlet.has_data_from(item.lower()):
                                try:
                                    shell.process_data()
                                    inlet.add_data_from_shell(shell)
                                except Exception:
                                    logging.exception(
                                        f"Error encountered when processing {file_name}"
                                    )
                                    continue
            if "HISTORY" in dirs:
                dirs.remove("HISTORY")

        # hakai data
        if from_csv:
            for file in fnmatch.filter(os.listdir(data_dir), "*.csv"):
                data = pandas.read_csv(os.path.join(data_dir, file))
                for inlet in inlet_list:
                    inside_inlet = data.loc[
                        lambda frame: [
                            inlet.contains(longitude=lon, latitude=lat)
                            for lon, lat in zip(frame["Longitude"], frame["Latitude"])
                        ]
                    ]
                    if len(inside_inlet) == 0:
                        continue
                    inlet.add_data_from_csv(inside_inlet, file)

    return inlet_list


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    parser.add_argument("-e", "--from-erddap", action="store_true")
    parser.add_argument("-n", "--from-netcdf", action="store_true")
    parser.add_argument("-c", "--from-csv", action="store_true")
    args = parser.parse_args()
    print("Preparing database")
    get_inlets(
        args.data,
        from_saved=False,
        from_netcdf=args.from_netcdf,
        from_erddap=args.from_erddap,
        from_csv=args.from_csv,
    )


if __name__ == "__main__":
    main()
