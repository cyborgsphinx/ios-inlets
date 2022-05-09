from dataclasses import dataclass
import logging
from typing import List

import convert
from erddapy import ERDDAP
import gsw
import numpy
import pandas
import urllib

CHUNKSIZE = 100000

UNASSIGNED = 0
UNALTERED = 1
COMPUTED = 2
ASSUMED = 3

@dataclass
class VariableData:
    ioos_category: str
    units: str
    standard_names: List[str]
    long_names: List[str] # some datasets do not have properly configured standard names
    fallback_names: List[str] # sometimes the erddap server is poorly configured and running locally

VARIABLES = {
    "required": {
        "latitude": VariableData(
            ioos_category="latitude",
            units="degrees_north",
            standard_names=["latitude"],
            long_names=["latitude"],
            fallback_names=[],
        ),
        "longitude": VariableData(
            ioos_category="longitude",
            units="degrees_east",
            standard_names=["longitude"],
            long_names=["longitude"],
            fallback_names=[],
        ),
        "time": VariableData(
            ioos_category="time",
            units="UTC",
            standard_names=["time"],
            long_names=["time"],
            fallback_names=[],
        ),
        "depth": VariableData(
            ioos_category="depth",
            units="m",
            standard_names=["depth", "instrument_depth"],
            long_names=["depth", "instrument_depth"],
            fallback_names=[],
        ),
        "source": VariableData(
            ioos_category="identifier",
            units="",
            standard_names=["file_name"],
            long_names=["filename"],
            fallback_names=["filename"],
        )
    },
    "optional": {
        "temperature": VariableData(
            ioos_category="temperature",
            units="C",
            standard_names=[
                "sea_water_temperature",
                "sea_water_conservative_temperature",
                "sea_water_potential_temperature",
            ],
            long_names=[
                "adcp_transducer_temp.",
                "sea_water_temperature",
            ],
            fallback_names=[
                "TEMPRTN1",
                "TEMPST01",
                "TEMPPR01",
                "TEMPPR03",
                "TEMPS901",
                "TEMPS601",
            ],
        ),
        "salinity": VariableData(
            ioos_category="salinity",
            units="PSU",
            standard_names=[
                "sea_water_salinity",
                "sea_water_practical_salinity",
                "sea_water_absolue_salinity",
                "sea_water_reference_salinity",
            ],
            long_names=[
                "sea_water_practical_salinity",
            ],
            fallback_names=[
                "PSLTZZ01",
                "ODSDM021",
                "SSALST01",
                "PSALST01",
                "PSALBST1",
                "sea_water_practical_salinity",
            ],
        ),
        "oxygen": VariableData(
            ioos_category="dissolved_o2",
            units="mL/L",
            standard_names=[
                "dissolved_oxygen_concentration",
                "volume_fraction_of_oxygen_in_water",
                "fractional_saturation_of_oxygen_in_water",
                "mass_concentration_of_oxygen_in_sea_water",
                "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water",
                "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_at_saturation",
                "moles_of_oxygen_per_unit_mass_in_sea_water",
            ],
            long_names=[
                "oxygen_concentration",
            ],
            fallback_names=["DOXYZZ01", "DOXMZZ01"],
        ),
    },
    "extra": {
        "pressure": VariableData(
            ioos_category="pressure",
            units="dbar",
            standard_names=["sea_water_pressure"],
            long_names=["sea_water_pressure_in_dbar"],
            fallback_names=[],
        ),
        "density": VariableData(
            ioos_category="",
            units="any",
            standard_names=[
                "sea_water_density",
                "sea_water_potential_density",
                "sea_water_neutral_density",
                "sea_water_sigma_t",
                "sea_water_sigma_theta",
            ],
            long_names=[],
            fallback_names=[],
        ),
    }
}

SERVERS = [
    #"http://localhost:8080", # for testing
    "https://data.cioospacific.ca", # CIOOS pacific, including IOS
    "https://catalogue.hakai.org", # Hakai
    #"https://salishsea.eos.ubc.ca", # UBC
]

CELSIUS_SPELLINGS = [
    name.lower()
    for name in [
        "C",
        "degC",
        "degree C",
        "degrees C",
        "deg_C",
        "degree_C",
        "degrees_C",
    ]
]

PSU_SPELLINGS = [
    name.lower()
    for name in [
        "PSU",
        "PSS-78",
    ]
]

DBAR_SPELLINGS = [
    name.lower()
    for name in [
        "dbar",
        "decibar",
    ]
]

UMOL_KG_SPELLINGS = [
    name.lower()
    for name in [
        "umol/kg",
        "μmole/kg",
        "mol/m",
        "mol/m**3",
        "\\u00ce\\u00bcmole/kg",
        "\u00ce\u00bcmole/kg",
        "\\u00c2\\u00b5mole/kg",
        "\u00c2\u00b5mole/kg",
    ]
]

def standardize_units(units):
    name = units.lower()
    return (
        "C" if name in CELSIUS_SPELLINGS
        else "PSU" if name in PSU_SPELLINGS
        else "dbar" if name in DBAR_SPELLINGS
        else "umol/kg" if name in UMOL_KG_SPELLINGS
        else units
    )


def search_to_download_key(key):
    # match statement before 3.10
    try:
        return {
            "min_lon": "longitude>=",
            "max_lon": "longitude<=",
            "min_lat": "latitude>=",
            "max_lat": "latitude<=",
        }[key]
    except KeyError:
        # default
        return key


def search_to_download(table):
    return {
        search_to_download_key(k): v
        for k, v in table.items()
    }


def convert_salinity(salinity, units):
    if units.lower() in ["psu", "pss-78"]:
        return salinity, UNALTERED
    elif units.lower() in ["ppt"]:
        return gsw.SP_from_SK(salinity), COMPUTED
    elif units.lower() in ["umol/kg"]:
        g_per_umol = 58.44 / 1000 / 1000
        return gsw.SP_from_SR(salinity * g_per_umol), COMPUTED
    else:
        logging.warning(f"Unknown units: {units}")
        return salinity, UNASSIGNED


def convert_oxygen(oxygen, units, data):
    if units.lower() in ["ml/l"]:
        return oxygen, UNALTERED
    elif units.lower() in ["mg/l"]:
        oxygen_mg_per_ml = 1.429
        return oxygen / oxygen_mg_per_ml, COMPUTED
    elif units.lower() in ["umol/kg", "mmol/m**3", "μmole/kg", "\\u00ce\\u00bcmole/kg",  "\\u00c2\\u00b5mole/kg"]:
        out, assumed = convert.convert_umol_kg_to_mL_L(
            oxygen,
            data["longitude"],
            data["latitude"],
            data["aggregated_temperature"],
            data["aggregated_salinity"],
            data["aggregated_pressure"],
        )
        return out, ASSUMED if assumed else COMPUTED
    elif units.lower() in ["umol/l", "μmole/l", "\\u00ce\\u00bcmole/l"]:
        oxygen_umol_per_ml = 44.661
        return oxygen / oxygen_umol_per_ml, COMPUTED
    else:
        logging.warning(f"Unknown units: {units}")
        return oxygen, UNASSIGNED


def combine_columns(data, desired_unit, units, new_column, columns, convert_fn, default=numpy.nan):
    new_metadata = new_column + "_metadata"
    new_quality = new_column + "_quality"
    new_data = data.assign(**{
        new_column: lambda _: default,
        new_metadata: lambda _: UNASSIGNED,
        new_quality: lambda _: 0,
    })
    standard_desired_unit = standardize_units(desired_unit)
    for column in columns:
        if len(desired_unit) == 0 or standardize_units(units.at[0, column]) == standard_desired_unit:
            # overwrite any converted data with something trusted
            if isinstance(default, float):
                indexer = new_data[new_column].isna() & (new_data[new_metadata] != UNALTERED)
            else:
                indexer = new_data[new_column] == default
            new_data.loc[indexer, new_column] = new_data[column]
            new_data.loc[indexer, new_metadata] = UNALTERED
        else:
            if isinstance(default, float):
                indexer = new_data[new_column].isna() & data[column].notna()
            else:
                indexer = (new_data[new_column] == default) & (data[column] != default)
            # expect convert_fn to return a (data, metadata) pair
            column_data, column_metadata = convert_fn(
                new_data.loc[indexer, column],
                units.at[0, column],
                new_data.loc[indexer],
            )
            new_data.loc[indexer, new_column] = column_data
            new_data.loc[indexer, new_metadata] = column_metadata
    return new_data


def process_data(data, variable_map, units, dataset_name):
    with_source = combine_columns(
        data,
        VARIABLES["required"]["source"].units,
        units,
        "source",
        variable_map["source"],
        lambda x, _units, _df: (x, UNALTERED),
        default=dataset_name,
    )
    with_pressure = combine_columns(
        with_source,
        VARIABLES["extra"]["pressure"].units,
        units,
        "aggregated_pressure",
        variable_map["pressure"],
        lambda x, _units, _df: (x, UNALTERED),
    )
    with_temp = combine_columns(
        with_pressure,
        VARIABLES["optional"]["temperature"].units,
        units,
        "aggregated_temperature",
        variable_map["temperature"],
        lambda x, _units, _df: (x, UNALTERED),
    )
    with_salinity = combine_columns(
        with_temp,
        VARIABLES["optional"]["salinity"].units,
        units,
        "aggregated_salinity",
        variable_map["salinity"],
        lambda x, from_units, df: convert_salinity(x, from_units)
    )
    with_oxygen = combine_columns(
        with_salinity,
        VARIABLES["optional"]["oxygen"].units,
        units,
        "aggregated_oxygen",
        variable_map["oxygen"],
        lambda x, from_units, df: convert_oxygen(x, from_units, df)
    )
    return with_oxygen


def find_variables_for(variable_data, server, dataset):
    found_names = server.get_var_by_attr(dataset_id=dataset, ioos_category=variable_data.ioos_category)
    for standard_name in variable_data.standard_names:
        found_names += server.get_var_by_attr(dataset_id=dataset, standard_name=standard_name)
    for long_name in variable_data.long_names:
        found_names += server.get_var_by_attr(dataset_id=dataset, long_name=long_name)
    if len(found_names) == 0:
        info = pandas.read_csv(server.get_info_url(dataset_id=dataset, response="csv"))
        found_names += [
            name
            for name in info.loc[info["Row Type"] == "variable", "Variable Name"]
            if name in variable_data.fallback_names
        ]
    return list(set(found_names))


def pull_data_for(inlet):
    for server in SERVERS:
        server_url = f"{server}/erddap"
        logging.info("Reading data from", server_url)
        e = ERDDAP(server=server_url, protocol="tabledap")
        parameters = inlet.bounding_box()
        search_url = e.get_search_url(response="csv", **parameters)
        try:
            results = pandas.read_csv(search_url)
        except urllib.error.HTTPError as err:
            if err.code in [500]:
                logging.warning(f"{server_url} is unavailable for searching")
                continue
            else:
                logging.error(f"{search_url} returned error {err.code}: {err.reason}")
                raise err
        for dataset in results["Dataset ID"].values:
            info_url = e.get_info_url(dataset_id=dataset, response="csv")
            try:
                info = pandas.read_csv(info_url)
            except urllib.error.HTTPError as err:
                if err.code in [500]:
                    logging.warning(f"{server_url} failed to retrieve info for {dataset}")
                    continue
                else:
                    logging.error(f"{info_url} returned error {err.code}: {err.reason}")
                    raise err
            logging.info(f"reading from {server_url}: {dataset}")
            name_map = {}

            for name, variable in VARIABLES["required"].items():
                name_map[name] = find_variables_for(variable, e, dataset)
                logging.debug(f"{name}: {name_map[name]}")
            if not all(
                len(name_map[name]) > 0
                for name in VARIABLES["required"].keys()
            ):
                logging.info(f"Missing critical information for {inlet.name} in {dataset} {VARIABLES['required'].keys()}")
                continue

            for name, variable in VARIABLES["optional"].items():
                name_map[name] = find_variables_for(variable, e, dataset)
                logging.debug(f"{name}: {name_map[name]}")
            if not any(
                len(name_map[name]) > 0
                for name in VARIABLES["optional"].keys()
            ):
                logging.info(f"No relevant data for {inlet.name} in {dataset} {VARIABLES['optional'].keys()}")
                continue

            for name, variable in VARIABLES["extra"].items():
                name_map[name] = find_variables_for(variable, e, dataset)
                logging.debug(f"{name}: {name_map[name]}")

            usable_variables = [item for value in name_map.values() for item in value]
            dl = e.get_download_url(
                dataset_id=dataset,
                variables=usable_variables,
                response="csv",
                constraints=search_to_download(parameters),
            )
            base, query = dl.split("?")
            dl = "?".join((base, urllib.parse.quote(query)))
            try:
                units = pandas.read_csv(dl, nrows=1)
                with pandas.read_csv(dl, chunksize=CHUNKSIZE, skiprows=(1,)) as reader:
                    for chunk in reader:
                        yield process_data(chunk, name_map, units, dataset)
            except urllib.error.HTTPError as err:
                if err.code in [404]:
                    logging.info(f"{dataset} turned up in advanced search despite not having data from {inlet.name}")
                else:
                    logging.error(f"{dl} returned error {err.code}: {err.reason}")
                    raise err


if __name__ == "__main__":
    main()
