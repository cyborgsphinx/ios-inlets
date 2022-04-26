from dataclasses import dataclass
import datetime
import logging
import os
import sqlite3
from typing import List


sqlite3.paramstyle = "named"
DB_NAME = os.path.join("data", "inlet_data.db")
SURFACE = "surface"
SHALLOW = "shallow"
DEEP = "deep"
DEEPER = "deeper"
DEEPEST = "deepest"
IGNORE = "ignore"
ALL = "all"
USED = "used"
USED_SURFACE = "used surface"
USED_DEEP = "used deep"


def _table_name(inlet_name: str) -> str:
    return inlet_name.lower().replace(" ", "_")


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


def _averaged(data: List[InletData]) -> List[InletData]:
    """Perform daily vertical averaging inside depth categories."""
    freqs = {}
    for datum in data:
        filename = datum.filename
        date = datum.time.date()
        key = (filename, date)
        if key not in freqs:
            freqs[key] = (0, 0)
        total, count = freqs[key]
        freqs[key] = (total + datum.value, count + 1)
    return [
        InletData(d, USED, 0, total / count, 0, 0, 0, fn)
        for (fn, d), (total, count) in freqs.items()
    ]


class InletDb:
    def __init__(self, inlet_name: str, clear: bool = False, db_name: str = DB_NAME):
        self.name = _table_name(inlet_name)
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
            self.__add_value(value, "temperature")
        except sqlite3.IntegrityError:
            logging.exception(
                f"Integrity error inserting temperature data ({value}) into database for {self.name}"
            )

    def add_temperature_data(self, data: List[InletData]):
        try:
            self.__add_data(data, "temperature")
        except sqlite3.IntegrityError:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting temperature data from {filename} into database for {self.name}"
            )

    def get_temperature_data(self, bucket, average: bool = False) -> List[InletData]:
        data = self.__get_data("temperature", bucket)
        if average:
            return _averaged(data)
        else:
            return data

    def add_salinity_value(self, value: InletData):
        try:
            self.__add_value(value, "salinity")
        except sqlite3.IntegrityError:
            logging.exception(
                f"Integrity error inserting salinity data ({value}) into database for {self.name}"
            )

    def add_salinity_data(self, data: List[InletData]):
        try:
            self.__add_data(data, "salinity")
        except sqlite3.IntegrityError:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting salinity data from {filename} into database for {self.name}"
            )

    def get_salinity_data(self, bucket, average: bool = False) -> List[InletData]:
        data = self.__get_data("salinity", bucket)
        if average:
            return _averaged(data)
        else:
            return data

    def add_oxygen_value(self, value: InletData):
        try:
            self.__add_value(value, "oxygen")
        except sqlite3.IntegrityError:
            logging.exception(
                f"Integrity error inserting oxygen data ({value}) into database for {self.name}"
            )

    def add_oxygen_data(self, data: List[InletData]):
        try:
            self.__add_data(data, "oxygen")
        except sqlite3.IntegrityError:
            filename = data[0].filename
            logging.exception(
                f"Integrity error inserting oxygen data from {filename} into database for {self.name}"
            )

    def get_oxygen_data(self, bucket, average: bool = False) -> List[InletData]:
        data = self.__get_data("oxygen", bucket)
        if average:
            return _averaged(data)
        else:
            return data

    def __add_value(self, value: InletData, kind: str):
        with self.connection:
            self.connection.execute(
                f"""
                insert into {self.name}
                values (
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
                {"kind": kind, **value.as_dict()},
            )

    def __add_data(self, data: List[InletData], kind: str):
        with self.connection:
            self.connection.executemany(
                f"""
                insert into {self.name}
                values (
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
                ({"kind": kind, **datum.as_dict()} for datum in data),
            )

    def __get_data(self, kind: str, bucket: str) -> List[InletData]:
        if bucket == ALL:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind
                """,
                {"kind": kind},
            )
        elif bucket == USED:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and bucket!=:bucket
                """,
                {"kind": kind, "bucket": IGNORE}
            )
        elif bucket == USED_SURFACE:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and bucket in (:surface, :shallow)
                """,
                {"kind": kind, "surface": SURFACE, "shallow": SHALLOW}
            )
        elif bucket == USED_DEEP:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and bucket in (:deep, :deeper, :deepest)
                """,
                {"kind": kind, "deep": DEEP, "deeper": DEEPER, "deepest": DEEPEST}
            )
        else:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and bucket=:bucket
                """,
                {"kind": kind, "bucket": bucket},
            )
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
            for row in cursor
        ]

    def __ensure_data_table(self):
        if not self.__has_data_table():
            with self.connection:
                self.connection.execute(
                    f"""
                    create table {self.name} (
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
                        assumed_density integer not null
                    )"""
                )

    def __clear_data_table(self):
        if self.__has_data_table():
            with self.connection:
                self.connection.execute(f"""drop table {self.name}""")

    def __has_data_table(self):
        cursor = self.connection.execute(
            f"""
            select count(name)
            from sqlite_master
            where type='table' and name='{self.name}'
            """
        )
        return cursor.fetchone()[0] > 0
