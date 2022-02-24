from dataclasses import dataclass
import datetime
import logging
import os
import sqlite3
from typing import Any, Dict, List, Set


sqlite3.paramstyle = "named"
DB_NAME = os.path.join("data", "inlet_data.db")


@dataclass(frozen=True)
class InletData:
    time: datetime.datetime
    bucket: str
    depth: float
    value: float
    quality: int
    longitude: float
    latitude: float
    filename: str
    computed: bool = False
    assumed_density: bool = False

    def as_dict(self):
        return {
            "time": self.time.isoformat(timespec="microseconds"),
            "bucket": self.bucket,
            "depth": self.depth,
            "value": self.value,
            "quality": self.quality,
            "longitude": self.longitude,
            "latitude": self.latitude,
            "filename": self.filename,
            "computed": self.computed,
            "assumed_density": self.assumed_density,
        }


class InletDb:
    def __init__(self, inlet_name: str, clear: bool = False, db_name: str = DB_NAME):
        self.name = inlet_name
        self.connection = sqlite3.connect(db_name)
        self.connection.row_factory = sqlite3.Row
        if clear:
            self.__clear_data_table()
        self.__ensure_data_table()

    def __del__(self):
        self.connection.close()

    def clear(self):
        self.__clear_data_table()

    def add_temperature_value(self, value: InletData):
        try:
            self.__add_value({"kind": "temperature", **value.as_dict()})
        except sqlite3.IntegrityError as e:
            logging.exception(
                f"Integrity error inserting temperature data ({value}) into database for {self.name}:\n{e}"
            )

    def add_temperature_data(self, data: List[InletData]):
        try:
            self.__add_data(set(data), {"kind": "temperature"})
        except sqlite3.IntegrityError as e:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting temperature data from {filename} into database for {self.name}:\n{e}"
            )
            logging.exception(f"{data}")

    def get_temperature_data(self) -> List[InletData]:
        return self.__get_data("temperature")

    def add_salinity_value(self, value: InletData):
        try:
            self.__add_value({"kind": "salinity", **value.as_dict()})
        except sqlite3.IntegrityError as e:
            logging.exception(
                f"Integrity error inserting salinity data ({value}) into database for {self.name}:\n{e}"
            )

    def add_salinity_data(self, data: List[InletData]):
        try:
            self.__add_data(set(data), {"kind": "salinity"})
        except sqlite3.IntegrityError as e:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting salinity data from {filename} into database for {self.name}:\n{e}"
            )
            logging.exception(f"{data}")

    def get_salinity_data(self) -> List[InletData]:
        return self.__get_data("salinity")

    def add_oxygen_value(self, value: InletData):
        try:
            self.__add_value({"kind": "oxygen", **value.as_dict()})
        except sqlite3.IntegrityError as e:
            logging.exception(
                f"Integrity error inserting oxygen data ({value}) into database for {self.name}:\n{e}"
            )

    def add_oxygen_data(self, data: List[InletData]):
        try:
            self.__add_data(set(data), {"kind": "oxygen"})
        except sqlite3.IntegrityError as e:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting oxygen data from {filename} into database for {self.name}:\n{e}"
            )
            logging.exception(f"{data}")

    def get_oxygen_data(self) -> List[InletData]:
        return self.__get_data("oxygen")

    def __add_value(self, value: Dict[str, Any]):
        with self.connection:
            self.connection.execute(
                """
                insert into data
                values (
                    :name,
                    :kind,
                    :filename,
                    :latitude,
                    :longitude,
                    :time,
                    :depth,
                    :value,
                    :bucket,
                    :quality,
                    :computed,
                    :assumed_density
                )""",
                {"name": self.name, **value},
            )

    def __add_data(self, data: Set[InletData], dict: Dict[str, Any]):
        with self.connection:
            self.connection.executemany(
                """
                insert into data
                values (
                    :name,
                    :kind,
                    :filename,
                    :latitude,
                    :longitude,
                    :time,
                    :depth,
                    :value,
                    :bucket,
                    :quality,
                    :computed,
                    :assumed_density
                )""",
                ({"name": self.name, **dict, **datum.as_dict()} for datum in data),
            )

    def __get_data(self, kind: str) -> List[InletData]:
        return [
            InletData(
                filename=row["filename"],
                latitude=row["latitude"],
                longitude=row["longitude"],
                time=datetime.datetime.fromisoformat(row["time"]),
                depth=row["depth"],
                value=row["value"],
                bucket=row["bucket"],
                quality=row["quality"],
                computed=(row["computed"] > 0),
                assumed_density=(row["assumed_density"] > 0),
            )
            for row in self.connection.execute(
                """select * from data
                where name=:name and kind=:kind
                """,
                {"name": self.name, "kind": kind},
            )
        ]

    def __ensure_data_table(self):
        if not self.__has_data_table():
            with self.connection:
                self.connection.execute(
                    """create table data (
                        name text not null,
                        kind text not null,
                        filename text not null,
                        latitude real not null,
                        longitude real not null,
                        time text not null,
                        depth real not null,
                        value real not null,
                        bucket text not null,
                        quality integer not null,
                        computed integer not null,
                        assumed_density integer not null,
                        primary key (
                            name,
                            kind,
                            filename,
                            time,
                            depth
                        )
                    )"""
                )

    def __clear_data_table(self):
        if self.__has_data_table():
            with self.connection:
                self.connection.execute("""drop table data""")

    def __has_data_table(self):
        cursor = self.connection.execute(
            """select count(name)
            from sqlite_master
            where type='table' and name='data'
            """
        )
        return cursor.fetchone()[0] > 0
