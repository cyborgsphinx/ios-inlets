from dataclasses import dataclass
from typing import List

from erddapy import ERDDAP
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
    "http://localhost:8080", # for testing
    #"https://data.cioospacific.ca", # CIOOS pacific, including IOS
    #"https://catalogue.hakai.org", # Hakai
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

def standardize_units(units):
    name = units.lower()
    return (
        "C" if name in CELSIUS_SPELLINGS
        else "PSU" if name in PSU_SPELLINGS
        else "dbar" if name in DBAR_SPELLINGS
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
        print(f"Unknown units: {units}")
        return salinity, UNALTERED


def combine_columns(data, desired_unit, units, new_column, columns, convert_fn):
    new_metadata = new_column + "_metadata"
    nulls = numpy.full(len(data), numpy.nan)
    emptys = numpy.full(len(data), UNASSIGNED)
    new_data = data.assign(**{new_column: pandas.Series(nulls), new_metadata: pandas.Series(emptys)})
    standard_desired_unit = standardize_units(desired_unit)
    for column in columns:
        if standardize_units(units.at[0, column]) == standard_desired_unit:
            # overwrite any converted data with something trusted
            indexer = new_data[new_column].isna() & (new_data[new_metadata] != UNALTERED)
            new_data.loc[indexer, new_column] = new_data[column]
            new_data.loc[indexer, new_metadata] = UNALTERED
        else:
            indexer = new_data[new_column].isna() & data[column].notna()
            # expect convert_fn to return a (data, metadata) pair
            column_data, column_metadata = convert_fn(new_data[column], units.at[0, column], new_data)
            new_data.loc[indexer, new_column] = column_data
            new_data.loc[indexer, new_metadata] = column_metadata
    return new_data


def process_data(data, variable_map, units):
    with_temp = combine_columns(
        data,
        VARIABLES["optional"]["temperature"].units,
        units,
        "aggregated_temperature",
        variable_map["temperature"],
        lambda x, units, df: (x, UNALTERED),
    )
    with_salinity = combine_columns(
        with_temp,
        VARIABLES["optional"]["salinity"].units,
        units,
        "aggregated_salinity",
        variable_map["salinity"],
        lambda x, from_units, df: convert_salinity(x, from_units)
    )
    print(with_salinity)


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


def pull_data_for(inlet_list):
    for server in SERVERS:
        server_url = f"{server}/erddap"
        print("Reading data from", server_url)
        e = ERDDAP(server=server_url, protocol="tabledap")
        for inlet in inlet_list:
            parameters = inlet.bounding_box()
            search_url = e.get_search_url(response="csv", **parameters)
            results = pandas.read_csv(search_url)
            for dataset in results["Dataset ID"].values:
                info_url = e.get_info_url(dataset_id=dataset, response="csv")
                info = pandas.read_csv(info_url)
                print(dataset)
                name_map = {}

                for name, variable in VARIABLES["required"].items():
                    name_map[name] = find_variables_for(variable, e, dataset)
                    print(f"{name}: {name_map[name]}")
                if not all(
                    len(name_map[name]) > 0
                    for name in VARIABLES["required"].keys()
                ):
                    print(f"Missing critical information for {inlet.name} in {dataset} {VARIABLES['required'].keys()}")
                    continue

                for name, variable in VARIABLES["optional"].items():
                    name_map[name] = find_variables_for(variable, e, dataset)
                    print(f"{name}: {name_map[name]}")
                if not any(
                    len(name_map[name]) > 0
                    for name in VARIABLES["optional"].keys()
                ):
                    print(f"No relevant data for {inlet.name} in {dataset} {VARIABLES['optional'].keys()}")
                    continue

                for name, variable in VARIABLES["extra"].items():
                    name_map[name] = find_variables_for(variable, e, dataset)
                    print(f"{name}: {name_map[name]}")

                usable_variables = [item for value in name_map.values() for item in value]
                dl = e.get_download_url(
                    dataset_id=dataset,
                    variables=["filename"] + usable_variables,
                    response="csv",
                    constraints=search_to_download(parameters),
                )
                base, query = dl.split("?")
                dl = "?".join((base, urllib.parse.quote(query)))
                print(dl)
                try:
                    units = pandas.read_csv(dl, nrows=1)
                    print(units)
                    with pandas.read_csv(dl, chunksize=CHUNKSIZE, skiprows=(1,)) as reader:
                        for chunk in reader:
                            process_data(chunk, name_map, units)
                except urllib.error.HTTPError as err:
                    if err.code == 404:
                        print(f"{dataset} turned up in advanced search despite not having data from {inlet.name}")
                    else:
                        raise err


def main():
    import inlets
    inlet_list = inlets.get_inlets("data", from_saved=True)
    pull_data_for(inlet_list[:1])


if __name__ == "__main__":
    main()
