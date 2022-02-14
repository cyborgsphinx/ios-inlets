import csv
from dataclasses import dataclass
import datetime
import fnmatch
import gsw
import ios_shell as ios
import logging
import math
import numpy
import os
from shapely.geometry import Polygon
import sqlite3
import xarray

import inlets

DB_NAME = os.path.join("data", "station_data.db")
CREATE = """create table data (
    filename text not null,
    latitude real not null,
    longitude real not null,
    time text not null,
    depth real,
    temperature real,
    temperature_quality integer,
    temperature_units text,
    salinity real,
    salinity_quality integer,
    salinity_units text,
    oxygen real,
    oxygen_quality integer,
    oxygen_units text,
    pressure real,
    pressure_quality integer,
    pressure_units text,
    primary key (
        filename,
        latitude,
        longitude,
        time,
        depth
    )
)"""

sqlite3.paramstyle = "named"


@dataclass
class StationData:
    filename: str
    latitude: float
    longitude: float
    time: datetime.datetime
    depth: float
    temperature_data: float = math.nan
    temperature_quality: int = -1
    temperature_units: str = ""
    salinity_data: float = math.nan
    salinity_quality: int = -1
    salinity_units: str = ""
    oxygen_data: float = math.nan
    oxygen_quality: int = -1
    oxygen_units: str = ""
    pressure_data: float = math.nan
    pressure_quality: int = -1
    pressure_units: str = ""

    def id(self):
        return self.filename + self.time.strftime("%Y/%m/%dT%H:%M:%S") + str(self.depth)

    def join(self, other):
        if self.id() != other.id():
            return
        logging.warning(
            f"{self.filename} has been determined to be from the same data measurement as {other.filename}"
        )

        if math.isnan(self.temperature_data):
            self.temperature_data = other.temperature_data
            self.temperature_quality = other.temperature_quality
            self.temperature_units = other.temperature_units
        if math.isnan(self.salinity_data):
            self.salinity_data = other.salinity_data
            self.salinity_quality = other.salinity_quality
            self.salinity_units = other.salinity_units
        if math.isnan(self.oxygen_data):
            self.oxygen_data = other.oxygen_data
            self.oxygen_quality = other.oxygen_quality
            self.oxygen_units = other.oxygen_units
        if math.isnan(self.pressure_data):
            self.pressure_data = other.pressure_data
            self.pressure_quality = other.pressure_quality
            self.pressure_units = other.pressure_units

    def as_dict(self):
        return {
            "filename": self.filename.split(".")[0],
            "latitude": self.latitude,
            "longitude": self.longitude,
            "time": self.time,
            "depth": self.depth,
            "temperature_data": self.temperature_data,
            "temperature_quality": self.temperature_quality,
            "temperature_units": self.temperature_units,
            "salinity_data": self.salinity_data,
            "salinity_quality": self.salinity_quality,
            "salinity_units": self.salinity_units,
            "oxygen_data": self.oxygen_data,
            "oxygen_quality": self.oxygen_quality,
            "oxygen_units": self.oxygen_units,
            "pressure_data": self.pressure_data,
            "pressure_quality": self.pressure_quality,
            "pressure_units": self.pressure_units,
        }


def _extend_data(arr, length):
    if arr is None:
        return numpy.full(length, numpy.nan)
    elif isinstance(arr, float):
        return numpy.full(length, arr)
    elif arr.size == 1:
        return numpy.full(length, arr.item())
    else:
        return inlets.get_array(arr)


def read_data(data_dir, inlet_list, skip_netcdf=False):
    good_files = set()
    if not skip_netcdf:
        for root, _, files in os.walk(os.path.join(data_dir, "netCDF_Data")):
            for item in fnmatch.filter(files, "*.nc"):
                print("NetCDF file:", item)
                file_name = os.path.join(root, item)
                data = xarray.open_dataset(file_name)
                good_file = True
                for inlet in inlet_list:
                    if inlet.contains(data):
                        assert any(dim in data.dims for dim in ["time", "z"])
                        if "time" in data.dims:
                            dim = inlets.get_array(data.time)
                        elif "z" in data.dims:
                            dim = inlets.get_array(data.z)
                        else:
                            raise ValueError(f"{item} has unknown dimensions")
                        length = dim.size
                        depth_data = _extend_data(inlets.find_depth_data(data), length)

                        temperature_data = inlets.find_temperature_data(data)
                        assert not isinstance(temperature_data, float)
                        temperature_units = (
                            temperature_data.units
                            if temperature_data is not None
                            else ""
                        )
                        temperature_data = _extend_data(temperature_data, length)

                        salinity_data = inlets.find_salinity_data(data)
                        assert not isinstance(salinity_data, float)
                        salinity_units = (
                            salinity_data.units if salinity_data is not None else ""
                        )
                        salinity_data = _extend_data(salinity_data, length)

                        oxygen_data = inlets.find_oxygen_data(data)
                        assert not isinstance(oxygen_data, float)
                        oxygen_units = (
                            oxygen_data.units if oxygen_data is not None else ""
                        )
                        oxygen_data = _extend_data(oxygen_data, length)

                        pressure_data = inlets.find_pressure_data(data)
                        assert not isinstance(pressure_data, float)
                        pressure_units = (
                            pressure_data.units if pressure_data is not None else ""
                        )
                        pressure_data = _extend_data(pressure_data, length)

                        for i in range(length):
                            if depth_data is not None:
                                depth_value = depth_data[i]
                            else:
                                good_file = False
                                depth_value = math.nan
                            if temperature_data is not None:
                                temperature_value = temperature_data[i]
                            else:
                                good_file = False
                                temperature_value = math.nan
                            if salinity_data is not None:
                                salinity_value = salinity_data[i]
                            else:
                                good_file = False
                                salinity_value = math.nan
                            if oxygen_data is not None:
                                oxygen_value = oxygen_data[i]
                            else:
                                good_file = False
                                oxygen_value = math.nan
                            if pressure_data is not None:
                                pressure_value = pressure_data[i]
                            else:
                                good_file = False
                                pressure_value = math.nan
                            yield StationData(
                                filename=item,
                                latitude=data.latitude.item(),
                                longitude=data.longitude.item(),
                                time=dim[i]
                                if "time" in data.dims
                                else data.time.item(),
                                depth=depth_value,
                                temperature_data=temperature_value,
                                temperature_units=temperature_units,
                                salinity_data=salinity_value,
                                salinity_units=salinity_units,
                                oxygen_data=oxygen_value,
                                oxygen_units=oxygen_units,
                                pressure_data=pressure_value,
                                pressure_units=pressure_units,
                            )
                if good_file:
                    good_files.add(item.split(".")[0])

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
                if item.split(".")[0] in good_files:
                    continue
                print("IOS Shell file:", item)
                good_file = True
                file_name = os.path.join(root, item)
                try:
                    shell = ios.ShellFile.fromfile(file_name, process_data=False)
                except Exception as e:
                    logging.exception(e)
                    continue
                for inlet in inlet_list:
                    if inlet.contains(shell.get_location()):
                        try:
                            shell.process_data()
                        except Exception as e:
                            logging.exception(e)
                            continue
                        channels = shell.file.channels
                        channel_details = shell.file.channel_details
                        names = [channel.name for channel in channels]

                        date_idx = inlets.find_column(channels, "Date")
                        time_idx = inlets.find_column(channels, "Time")

                        depth_idx = inlets.find_column(channels, "Depth", "m", "metre")
                        depth_pad = inlets.get_pad_value(channel_details, depth_idx)
                        if depth_pad is None or math.isnan(depth_pad):
                            depth_pad = -99

                        temperature_idx = inlets.find_column(
                            channels, "Temperature", "C", "'deg C'"
                        )
                        temperature_quality_idx = (
                            temperature_idx + 1
                            if inlets.has_quality(temperature_idx, names)
                            else -1
                        )
                        temperature_pad = inlets.get_pad_value(
                            channel_details, temperature_idx
                        )
                        if temperature_pad is None or math.isnan(temperature_pad):
                            temperature_pad = -99
                        temperature_units = inlets.get_units(channels, temperature_idx)

                        salinity_idx = inlets.find_column(
                            channels, "Salinity", "PSU", "PSS-78"
                        )
                        salinity_quality_idx = (
                            salinity_idx + 1
                            if inlets.has_quality(salinity_idx, names)
                            else -1
                        )
                        salinity_pad = inlets.get_pad_value(
                            channel_details, salinity_idx
                        )
                        if salinity_pad is None or math.isnan(salinity_pad):
                            salinity_pad = -99
                        salinity_units = inlets.get_units(channels, salinity_idx)

                        oxygen_idx = inlets.find_column(channels, "Oxygen", "mL/L")
                        oxygen_quality_idx = (
                            oxygen_idx + 1
                            if inlets.has_quality(oxygen_idx, names)
                            else -1
                        )
                        oxygen_pad = inlets.get_pad_value(channel_details, oxygen_idx)
                        if oxygen_pad is None or math.isnan(oxygen_pad):
                            oxygen_pad = -99
                        oxygen_units = inlets.get_units(channels, oxygen_idx)

                        pressure_idx = inlets.find_column(
                            channels, "Pressure", "dbar", "decibar"
                        )
                        pressure_quality_idx = (
                            pressure_idx + 1
                            if inlets.has_quality(pressure_idx, names)
                            else -1
                        )
                        pressure_pad = inlets.get_pad_value(
                            channel_details, pressure_idx
                        )
                        if pressure_pad is None or math.isnan(pressure_pad):
                            pressure_pad = -99
                        pressure_units = inlets.get_units(channels, pressure_idx)

                        for row in shell.data:
                            if date_idx >= 0:
                                date = row[date_idx]
                            else:
                                date = shell.file.start_time
                            if time_idx >= 0:
                                time = row[time_idx]
                                if isinstance(date, datetime.date) and not isinstance(
                                    date, datetime.datetime
                                ):
                                    assert isinstance(time, datetime.time)
                                    date = datetime.datetime.combine(
                                        date, time, tzinfo=shell.get_time().tzinfo
                                    )
                            assert isinstance(date, datetime.datetime)

                            if depth_idx >= 0:
                                depth = row[depth_idx]
                            elif shell.instrument is not None and not math.isnan(
                                shell.instrument.depth
                            ):
                                depth = shell.instrument.depth
                            else:
                                good_file = False
                                depth = math.nan
                            assert isinstance(depth, float) or isinstance(depth, int)

                            if temperature_idx >= 0:
                                temperature = row[temperature_idx]
                            else:
                                good_file = False
                                temperature = math.nan
                            assert isinstance(temperature, float) or isinstance(
                                temperature, int
                            )
                            if temperature_quality_idx >= 0:
                                try:
                                    temperature_quality = int(
                                        row[temperature_quality_idx]
                                    )
                                except ValueError:
                                    temperature_quality = -1
                            else:
                                temperature_quality = -1

                            if salinity_idx >= 0:
                                salinity = row[salinity_idx]
                            else:
                                good_file = False
                                salinity = math.nan
                            assert isinstance(salinity, float) or isinstance(
                                salinity, int
                            )
                            if salinity_quality_idx >= 0:
                                try:
                                    salinity_quality = int(row[salinity_quality_idx])
                                except ValueError:
                                    salinity_quality = -1
                            else:
                                salinity_quality = -1

                            if oxygen_idx >= 0:
                                oxygen = row[oxygen_idx]
                            else:
                                good_file = False
                                oxygen = math.nan
                            assert isinstance(oxygen, float) or isinstance(oxygen, int)
                            if oxygen_quality_idx >= 0:
                                try:
                                    oxygen_quality = int(row[oxygen_quality_idx])
                                except ValueError:
                                    oxygen_quality = -1
                            else:
                                oxygen_quality = -1

                            if pressure_idx >= 0:
                                pressure = row[pressure_idx]
                            else:
                                good_file = False
                                pressure = math.nan
                            assert isinstance(pressure, float) or isinstance(
                                pressure, int
                            )
                            if pressure_quality_idx >= 0:
                                try:
                                    pressure_quality = int(row[pressure_quality_idx])
                                except ValueError:
                                    pressure_quality = -1
                            else:
                                pressure_quality = -1

                            yield StationData(
                                filename=item,
                                latitude=shell.location.latitude,
                                longitude=shell.location.longitude,
                                time=date,
                                depth=depth,
                                temperature_data=temperature,
                                temperature_quality=temperature_quality,
                                temperature_units=temperature_units,
                                salinity_data=salinity,
                                salinity_quality=salinity_quality,
                                salinity_units=salinity_units,
                                oxygen_data=oxygen,
                                oxygen_quality=oxygen_quality,
                                oxygen_units=oxygen_units,
                                pressure_data=pressure,
                                pressure_quality=pressure_quality,
                                pressure_units=pressure_units,
                            )
                if good_file:
                    good_files.add(item.split(".")[0])

            if "HISTORY" in dirs:
                dirs.remove("HISTORY")

    # hakai data
    for file in fnmatch.filter(os.listdir(data_dir), "*.csv"):
        print("CSV file:", file)
        with open(os.path.join(data_dir, file)) as f:
            reader = csv.DictReader(f)
            for row in reader:
                coords = {
                    "longitude": float(row["Longitude"]),
                    "latitude": float(row["Latitude"]),
                }
                for inlet in inlet_list:
                    if inlet.contains(coords):
                        yield StationData(
                            filename=item,
                            latitude=coords["latitude"],
                            longitude=coords["longitude"],
                            time=datetime.datetime.fromisoformat(
                                row["Measurement time"]
                            ),
                            depth=float(row["Depth (m)"]),
                            temperature_data=float(row["Temperature (deg C)"]),
                            temperature_units="C",
                            salinity_data=float(row["Salinity (PSU)"]),
                            salinity_units="PSU",
                            oxygen_data=float(row["Dissolved O2 (mL/L)"]),
                            oxygen_units="mL/L",
                        )


class StationDb:
    def __init__(self, name):
        self.name = name
        self.connection = sqlite3.connect(name)
        self.cursor = self.connection.cursor()
        self.cursor.execute(
            """
            select count(name)
            from sqlite_master
            where type='table' and name='data'
        """
        )
        if self.cursor.fetchone()[0] == 0:
            self.cursor.execute(CREATE)

    def __del__(self):
        self.cursor.close()
        self.connection.close()

    def reset_data(self):
        self.cursor.execute("drop table data")
        self.cursor.execute(CREATE)

    def add_data(self, value: StationData):
        try:
            self.cursor.execute(
                """
            insert into data
            values (
                :filename,
                :latitude,
                :longitude,
                :time,
                :depth,
                :temperature_data,
                :temperature_quality,
                :temperature_units,
                :salinity_data,
                :salinity_quality,
                :salinity_units,
                :oxygen_data,
                :oxygen_quality,
                :oxygen_units,
                :pressure_data,
                :pressure_quality,
                :pressure_units
            )""",
                value.as_dict(),
            )
        except sqlite3.DatabaseError as e:
            self.cursor.execute(
                """
                select temperature, salinity, oxygen, pressure
                from data
                where filename=:filename and time=:time and depth=:depth
                """,
                value.as_dict(),
            )
            if (row := self.cursor.fetchone()) is not None:
                value_dict = value.as_dict()
                if row[0] is None or math.isnan(row[0]):
                    value_dict["value"] = value.temperature_data
                    value_dict["units"] = value.temperature_units
                    self.cursor.execute(
                        """
                        update data
                        set temperature = :value, temperature_units = :units
                        where filename=:filename and time=:time and depth=:depth
                    """,
                        value_dict,
                    )
                if row[1] is None or math.isnan(row[1]):
                    value_dict["value"] = value.salinity_data
                    value_dict["units"] = value.salinity_units
                    self.cursor.execute(
                        """
                        update data
                        set salinity = :value, salinity_units = :units
                        where filename=:filename and time=:time and depth=:depth
                    """,
                        value_dict,
                    )
                if row[2] is None or math.isnan(row[2]):
                    value_dict["value"] = value.oxygen_data
                    value_dict["units"] = value.oxygen_units
                    self.cursor.execute(
                        """
                        update data
                        set oxygen = :value, oxygen_units = :units
                        where filename=:filename and time=:time and depth=:depth
                    """,
                        value_dict,
                    )
                if row[3] is None or math.isnan(row[3]):
                    value_dict["value"] = value.pressure_data
                    value_dict["units"] = value.pressure_units
                    self.cursor.execute(
                        """
                        update data
                        set pressure = :value, pressure_units = :units
                        where filename=:filename and time=:time and depth=:depth
                    """,
                        value_dict,
                    )
            else:
                print(f"Unknown error: {e}")
        finally:
            self.connection.commit()

    def data(self):
        # iterator/generator over data in table
        self.cursor.execute("select * from data")
        data = self.cursor.fetchone()
        print(data)
        while data:
            yield data
            data = self.cursor.fetchone()


if __name__ == "__main__":
    import json

    inlet_list = []
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

    db = StationDb(DB_NAME)
    db.reset_data()

    for data in read_data("data", inlet_list, skip_netcdf=False):
        db.add_data(data)
    print("Done")
