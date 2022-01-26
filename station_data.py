import csv
from dataclasses import dataclass
import datetime
import fnmatch
import ios_shell as ios
import logging
import math
import os
from shapely.geometry import Polygon
import sqlite3
import xarray

import inlets

DB_NAME = os.path.join("data", "station_data.db")
CREATE = """create table data (
filename text,
latitude real,
longitude real,
kind text,
project text,
mission text,
station test,
event int,
time text,
depth real,
temperature real,
temperature_units text,
salinity real,
salinity_units text,
oxygen real,
oxygen_units text,
pressure real,
pressure_units text
)"""

sqlite3.paramstyle = "named"


@dataclass
class StationData():
    filename: str
    latitude: float
    longitude: float
    kind: str
    project: str
    mission: str
    station: str
    event: int
    time: datetime.datetime
    depth: float
    temperature_data: float = math.nan
    temperature_units: str = ""
    salinity_data: float = math.nan
    salinity_units: str = ""
    oxygen_data: float = math.nan
    oxygen_units: str = ""
    pressure_data: float = math.nan
    pressure_units: str = ""

    def id(self):
        return (str(self.latitude)
            + str(self.longitude)
            + str(self.kind)
            + self.project
            + self.mission
            + self.station
            + str(self.event)
            + self.time.strftime("%Y/%m/%dT%H:%M:%S")
            + str(self.depth))

    def join(self, other):
        if self.id() != other.id():
            return
        logging.warning(f"{self.filename} has been determined to be from the same data measurement as {other.filename}")

        if math.isnan(self.temperature_data):
            self.temperature_data = other.temperature_data
            self.temperature_units = other.temperature_units
        if math.isnan(self.salinity_data):
            self.salinity_data = other.salinity_data
            self.salinity_units = other.salinity_units
        if math.isnan(self.oxygen_data):
            self.oxygen_data = other.oxygen_data
            self.oxygen_units = other.oxygen_units
        if math.isnan(self.pressure_data):
            self.pressure_data = other.pressure_data
            self.pressure_units = other.pressure_units

    def as_dict(self):
        return {
            "filename": self.filename,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "kind": self.kind,
            "project": self.project,
            "mission": self.mission,
            "station": self.station,
            "event": self.event,
            "time": self.time,
            "depth": self.depth,
            "temperature_data": self.temperature_data,
            "temperature_units": self.temperature_units,
            "salinity_data": self.salinity_data,
            "salinity_units": self.salinity_units,
            "oxygen_data": self.oxygen_data,
            "oxygen_units": self.oxygen_units,
            "pressure_data": self.pressure_data,
            "pressure_units": self.pressure_units,
        }


def _extend_data(arr, length):
    if arr is None:
        return [math.nan] * length
    elif isinstance(arr, float):
        return [arr] * length
    else:
        return inlets.get_array(arr)


def _data_kind(data):
    filename = data.filename.item()
    if "instrumentKind" in data:
        return data.instrumentKind.lower()
    elif filename.endswith(".bot") or filename.endswith(".che"):
        return "bottle"
    elif filename.endswith(".adcp") or filename.endswith(".L1"):
        return "adcp"
    else:
        raise ValueError("Unknown kind: " + data.filename)


def read_data(data_dir, inlet_list, skip_netcdf=False):
    if not skip_netcdf:
        for root, _, files in os.walk(os.path.join(data_dir, "netCDF_Data")):
            for item in fnmatch.filter(files, "*.nc"):
                try:
                    file_name = os.path.join(root, item)
                    data = xarray.open_dataset(file_name)
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
                            kind = _data_kind(data)

                            temperature_data = inlets.find_temperature_data(data)
                            assert not isinstance(temperature_data, float)
                            temperature_units = temperature_data.units if temperature_data is not None else ""
                            temperature_data = _extend_data(temperature_data, length)

                            salinity_data = inlets.find_salinity_data(data)
                            assert not isinstance(salinity_data, float)
                            salinity_units = salinity_data.units if salinity_data is not None else ""
                            salinity_data = _extend_data(salinity_data, length)

                            oxygen_data = inlets.find_oxygen_data(data)
                            assert not isinstance(oxygen_data, float)
                            oxygen_units = oxygen_data.units if oxygen_data is not None else ""
                            oxygen_data = _extend_data(oxygen_data, length)

                            pressure_data = inlets.find_pressure_data(data)
                            assert not isinstance(pressure_data, float)
                            pressure_units = pressure_data.units if pressure_data is not None else ""
                            pressure_data = _extend_data(pressure_data, length)

                            for i in range(length):
                                depth_value = depth_data[i]
                                if temperature_data is not None:
                                    temperature_value = temperature_data[i]
                                else:
                                    temperature_value = math.nan
                                if salinity_data is not None:
                                    salinity_value = salinity_data[i]
                                else:
                                    salinity_value = math.nan
                                if oxygen_data is not None:
                                    oxygen_value = oxygen_data[i]
                                else:
                                    oxygen_value = math.nan
                                if pressure_data is not None:
                                    pressure_value = pressure_data[i]
                                else:
                                    pressure_value = math.nan
                                yield StationData(
                                    filename=item,
                                    latitude=data.latitude.item(),
                                    longitude=data.longitude.item(),
                                    kind=kind,
                                    project=data.project,
                                    mission=data.deployment_cruise_number,
                                    station=data.station,
                                    event=data.deployment_number,
                                    time=dim[i] if "time" in data.dims else data.time.item(),
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
                except Exception as e:
                    logging.warning(f"Exception in {item}: {e}")
                    raise e

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
                    if inlet.contains(data):
                        channels = shell.file.channels
                        channel_details = shell.file.channel_details

                        date_idx = inlets.find_column(shell.file.channels, "Date")
                        time_idx = inlets.find_column(channels, "Time")

                        depth_idx = inlets.find_column(channels, "Depth", "m", "metre")
                        depth_pad = inlets.get_pad_value(channel_details, depth_idx)
                        if depth_pad is None or math.isnan(depth_pad):
                            depth_pad = -99

                        temperature_idx = inlets.find_column(channels, "Temperature", "C", "'deg C'")
                        temperature_pad = inlets.get_pad_value(channel_details, temperature_idx)
                        if temperature_pad is None or math.isnan(temperature_pad):
                            temperature_pad = -99
                        temperature_units = inlets.get_units(channels, temperature_idx)

                        salinity_idx = inlets.find_column(channels, "Salinity", "PSU", "PSS-78")
                        salinity_pad = inlets.get_pad_value(channel_details, salinity_idx)
                        if salinity_pad is None or math.isnan(salinity_pad):
                            salinity_pad = -99
                        salinity_units = inlets.get_units(channels, salinity_idx)

                        oxygen_idx = inlets.find_column(channels, "Oxygen", "mL/L")
                        oxygen_pad = inlets.get_pad_value(channel_details, oxygen_idx)
                        if oxygen_pad is None or math.isnan(oxygen_pad):
                            oxygen_pad = -99
                        oxygen_units = inlets.get_units(channels, oxygen_idx)

                        pressure_idx = inlets.find_column(channels, "Pressure", "dbar", "decibar")
                        pressure_pad = inlets.get_pad_value(channel_details, pressure_idx)
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
                                if isinstance(date, datetime.date):
                                    assert isinstance(time, datetime.time)
                                    date = datetime.datetime.combine(date, time, tzinfo=shell.get_time().tzinfo)
                            assert isinstance(date, datetime.datetime)

                            if depth_idx >= 0:
                                depth = row[depth_idx]
                            elif shell.instrument is not None and not math.isnan(shell.instrument.depth):
                                depth = shell.instrument.depth
                            else:
                                depth = math.nan
                            assert isinstance(depth, float)

                            if temperature_idx >= 0:
                                temperature = row[temperature_idx]
                            else:
                                temperature = math.nan
                            assert isinstance(temperature, float)

                            if salinity_idx >= 0:
                                salinity = row[salinity_idx]
                            else:
                                salinity = math.nan
                            assert isinstance(salinity, float)

                            if oxygen_idx >= 0:
                                oxygen = row[oxygen_idx]
                            else:
                                oxygen = math.nan
                            assert isinstance(oxygen, float)

                            if pressure_idx >= 0:
                                pressure = row[pressure_idx]
                            else:
                                pressure = math.nan
                            assert isinstance(pressure, float)

                            yield StationData(
                                filename=item,
                                latitude=shell.location.latitude,
                                longitude=shell.location.longitude,
                                kind=shell.file.data_description.lower(),
                                project=shell.administration.project,
                                mission=shell.administration.mission,
                                station=shell.location.station,
                                event=shell.location.event_number,
                                time=date,
                                depth=depth,
                                temperature_data=temperature,
                                temperature_units=temperature_units,
                                salinity_data=salinity,
                                salinity_units=salinity_units,
                                oxygen_data=oxygen,
                                oxygen_units=oxygen_units,
                                pressure_data=pressure,
                                pressure_units=pressure_units,
                            )

            if "HISTORY" in dirs:
                dirs.remove("HISTORY")

        # hakai data
        for file in fnmatch.filter(os.listdir(data_dir), "*.csv"):
            with open(os.path.join(data_dir, file)) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    coords = {
                        "longitude": float(row["Longitude"]),
                        "latitude": float(row["Latitude"]),
                    }
                    for inlet in inlet_list:
                        if inlet.contains(data):
                            yield StationData(
                                filename=item,
                                latitude=coords["latitude"],
                                longitude=coords["longitude"],
                                kind="n/a",
                                project="Hakai",
                                mission=row["Cruise"],
                                station=row["Station"],
                                event=-1,
                                time=datetime.datetime.fromisoformat(row["Measurement time"]),
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
        self.cursor.execute("""
            select count(name)
            from sqlite_master
            where type='table' and name='data'
        """)
        if self.cursor.fetchone()[0] == 0:
            self.cursor.execute(CREATE)

    def __del__(self):
        self.cursor.close()
        self.connection.close()

    def reset_data(self):
        self.cursor.execute("drop table data")
        self.cursor.execute(CREATE)

    def add_data(self, value: StationData):
        self.cursor.execute("""
            select temperature, salinity, oxygen, pressure
            from data
            where latitude=:latitude and longitude=:longitude and kind=:kind and project=:project and station=:station and event=:event and time=:time and depth=:depth
            """, value.as_dict())
        if (row := self.cursor.fetchone()) is not None:
            # already in database
            value_dict = value.as_dict()
            if row[0] is None or math.isnan(row[0]):
                value_dict["value"] = value.temperature_data
                value_dict["units"] = value.temperature_units
                self.cursor.execute("""
                    update data
                    set temperature = :value, temperature_units = :units
                    where latitude=:latitude and longitude=:longitude and kind=:kind and project=:project and station=:station and event=:event and time=:time and depth=:depth
                """, value_dict)
            if row[1] is None or math.isnan(row[1]):
                value_dict["name"] = "salinity"
                value_dict["name_units"] = "salinity_units"
                value_dict["value"] = value.salinity_data
                value_dict["units"] = value.salinity_units
                self.cursor.execute("""
                    update data
                    set salinity = :value, salinity_units = :units
                    where latitude=:latitude and longitude=:longitude and kind=:kind and project=:project and station=:station and event=:event and time=:time and depth=:depth
                """, value_dict)
            if row[2] is None or math.isnan(row[2]):
                value_dict["name"] = "oxygen"
                value_dict["name_units"] = "oxygen_units"
                value_dict["value"] = value.oxygen_data
                value_dict["units"] = value.oxygen_units
                self.cursor.execute("""
                    update data
                    set oxygen = :value, oxygen_units = :units
                    where latitude=:latitude and longitude=:longitude and kind=:kind and project=:project and station=:station and event=:event and time=:time and depth=:depth
                """, value_dict)
            if row[3] is None or math.isnan(row[3]):
                value_dict["name"] = "pressure"
                value_dict["name_units"] = "pressure_units"
                value_dict["value"] = value.temperature_data
                value_dict["units"] = value.temperature_units
                self.cursor.execute("""
                    update data
                    set pressure = :value, pressure_units = :units
                    where latitude=:latitude and longitude=:longitude and kind=:kind and project=:project and station=:station and event=:event and time=:time and depth=:depth
                """, value_dict)
        else:
            self.cursor.execute("""
            insert into data
            values (
                :filename,
                :latitude,
                :longitude,
                :kind,
                :project,
                :mission,
                :station,
                :event,
                :time,
                :depth,
                :temperature_data,
                :temperature_units,
                :salinity_data,
                :salinity_units,
                :oxygen_data,
                :oxygen_units,
                :pressure_data,
                :pressure_units
            )""", value.as_dict())

    def data(self):
        # iterator/generator over data in table
        pass


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
