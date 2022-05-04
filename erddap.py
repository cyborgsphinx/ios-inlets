from erddapy import ERDDAP
import pandas
import urllib

CHUNKSIZE = 100000

VARIABLES = {
    "latitude": ["latitude"],
    "longitude": ["longitude"],
    "time": ["time"],
    "depth": ["depth"],
    "temperature": [
        "sea_water_temperature",
        "sea_water_conservative_temperature",
        "sea_water_potential_temperature",
    ],
    "salinity": [
        "sea_water_salinity",
        "sea_water_practical_salinity",
        "sea_water_absolue_salinity",
        "sea_water_reference_salinity",
    ],
    "oxygen": [
        "volume_fraction_of_oxygen_in_water",
        "fractional_saturation_of_oxygen_in_water",
        "mass_concentration_of_oxygen_in_sea_water",
        "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water",
        "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_at_saturation",
        "moles_of_oxygen_per_unit_mass_in_sea_water",
    ],
    "pressure": ["sea_water_pressure"],
    "density": [
        "sea_water_density",
        "sea_water_potential_density",
        "sea_water_neutral_density",
        "sea_water_sigma_t",
        "sea_water_sigma_theta",
    ],
}

SERVERS = [
    "data.cioospacific.ca", # CIOOS pacific, including IOS
    #"catalogue.hakai.org", # Hakai
    #"salishsea.eos.ubc.ca", # UBC
]

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


def combine_columns(data, units, new_column, columns):
    new_units = new_column + "_units"
    nulls = numpy.full(len(data), numpy.nan)
    emptys = numpy.full(len(data), "")
    data.assign(**{new_column: pandas.Series(nulls), new_units: pandas.Series(emptys)})
    for column in columns:
        indexer = data[new_column].isna() & data[column].notna()
        data.loc[indexer, new_column] = data[column]
        data.loc[indexer, new_units] = units[column]

def process_data(data):
    print(data)

def pull_data_for(inlet_list):
    for server in SERVERS:
        server_url = f"https://{server}/erddap"
        print("Reading data from", server_url)
        e = ERDDAP(server=server_url, protocol="tabledap")
        for inlet in inlet_list:
            parameters = inlet.bounding_box()
            search_url = e.get_search_url(response="csv", **parameters)
            print(search_url)
            results = pandas.read_csv(search_url)
            for dataset in results["Dataset ID"].values:
                print(dataset)
                usable_variables = []
                for name, standard_names in VARIABLES.items():
                    print("finding standard names for", name)
                    for standard_name in standard_names:
                        names = e.get_var_by_attr(dataset_id=dataset, standard_name=standard_name)
                        print(names)
                        usable_variables += names
                dl = e.get_download_url(
                    dataset_id=dataset,
                    variables=usable_variables,
                    response="csv",
                    constraints=search_to_download(parameters),
                )
                print(dl)
                try:
                    with pandas.read_csv(dl, chunksize=CHUNKSIZE) as reader:
                        for chunk in reader:
                            process_data(chunk)
                except urllib.error.HTTPError as err:
                    if err.code == 404:
                        print(f"{dataset} turned up in advanced search despite not having data from {inlet.name}")
                    else:
                        raise err

def main():
    import inlets
    inlet_list = inlets.get_inlets("data", from_saved=True)
    pull_data_for(inlet_list)

if __name__ == "__main__":
    main()
