from dataclasses import dataclass
import datetime
import logging
import os
import sqlite3
from typing import List, Tuple


sqlite3.paramstyle = "named"
DB_NAME = os.path.join("data", "inlet_data.db")


def _table_name(inlet_name: str) -> str:
    return inlet_name.lower().replace(" ", "_")


@dataclass(frozen=True)
class InletData:
    time: datetime.datetime
    depth: float
    value: float
    quality: int
    longitude: float
    latitude: float
    source: str
    computed: bool = False
    assumed_density: bool = False

    def as_dict(self):
        return {
            "time": self.time.isoformat(timespec="microseconds"),
            "depth": self.depth,
            "value": self.value,
            "quality": self.quality,
            "longitude": self.longitude,
            "latitude": self.latitude,
            "source": self.source,
            "computed": self.computed,
            "assumed_density": self.assumed_density,
        }


def _averaged(data: List[InletData]) -> List[InletData]:
    """Perform daily vertical averaging inside depth categories."""
    freqs = {}
    for datum in data:
        source = datum.source
        date = datum.time.date()
        key = (source, date)
        if key not in freqs:
            freqs[key] = (0, 0)
        total, count = freqs[key]
        freqs[key] = (total + datum.value, count + 1)
    return [
        InletData(d, 0, total / count, 0, 0, 0, fn)
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
            source = data[0].source
            logging.exception(
                f"Integrity error inserting temperature data from {source} into database for {self.name}"
            )

    def get_temperature_data(self, bucket: Tuple[float, float], average: bool = False) -> List[InletData]:
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
            source = data[0].source
            logging.exception(
                f"Integrity error inserting salinity data from {source} into database for {self.name}"
            )

    def get_salinity_data(self, bucket: Tuple[float, float], average: bool = False) -> List[InletData]:
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
            source = data[0].source
            logging.exception(
                f"Integrity error inserting oxygen data from {source} into database for {self.name}"
            )

    def get_oxygen_data(self, bucket: Tuple[float, float], average: bool = False) -> List[InletData]:
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
                    :source,
                    :latitude,
                    :longitude,
                    :time,
                    :depth,
                    :value,
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
                    :source,
                    :latitude,
                    :longitude,
                    :time,
                    :depth,
                    :value,
                    :quality,
                    :computed,
                    :assumed_density
                )""",
                ({"kind": kind, **datum.as_dict()} for datum in data),
            )

    def __get_data(self, kind: str, bucket: Tuple[float, float]) -> List[InletData]:
        min_depth, max_depth = bucket
        if min_depth is None and max_depth is None:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind
                """,
                {"kind": kind},
            )
        elif min_depth is None:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and depth<=:max
                """,
                {"kind": kind, "max": max_depth},
            )
        elif max_depth is None:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and depth>=:min
                """,
                {"kind": kind, "min": min_depth},
            )
        else:
            cursor = self.connection.execute(
                f"""select * from {self.name}
                where kind=:kind and depth>=:min and depth<=:max
                """,
                {"kind": kind, "min": min_depth, "max": max_depth},
            )
        return [
            InletData(
                source=row["source"],
                latitude=row["latitude"],
                longitude=row["longitude"],
                time=datetime.datetime.fromisoformat(row["time"]),
                depth=row["depth"],
                value=row["value"],
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
                        source text not null,
                        latitude real not null,
                        longitude real not null,
                        time text not null,
                        depth real not null,
                        value real not null,
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
