import argparse
import inlets


def take_before(inlet, before):
    inlet.temperature_data = [datum for datum in inlet.temperature_data if datum.time.year < before]
    inlet.salinity_data = [datum for datum in inlet.salinity_data if datum.time.year < before]
    inlet.oxygen_data = [datum for datum in inlet.oxygen_data if datum.time.year < before]
    return inlet


def take_greater_than(inlet, greater_than):
    inlet.temperature_data = [datum for datum in inlet.temperature_data if datum.datum > greater_than]
    inlet.salinity_data = [datum for datum in inlet.salinity_data if datum.datum > greater_than]
    inlet.oxygen_data = [datum for datum in inlet.oxygen_data if datum.datum > greater_than]
    return inlet


def take_lesser_than(inlet, lesser_than):
    inlet.temperature_data = [datum for datum in inlet.temperature_data if datum.datum < lesser_than]
    inlet.salinity_data = [datum for datum in inlet.salinity_data if datum.datum < lesser_than]
    inlet.oxygen_data = [datum for datum in inlet.oxygen_data if datum.datum < lesser_than]
    return inlet


def main():
    parser = argparse.ArgumentParser()
    # inlet retrieval args
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-n", "--skip-netcdf", action="store_true")
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    # search args
    parser.add_argument("-b", "--before", metavar="YEAR", type=int, nargs="?", default=None)
    parser.add_argument("-g", "--greater-than", type=int, nargs="?", default=None)
    parser.add_argument("-l", "--lesser-than", type=int, nargs="?", default=None)
    args = parser.parse_args()
    inlet_list = inlets.get_inlets(args.data, args.from_saved, args.skip_netcdf)

    if args.before is not None:
        inlet_list = [take_before(inlet, args.before) for inlet in inlet_list]
    if args.greater_than is not None:
        inlet_list = [take_greater_than(inlet, args.greater_than) for inlet in inlet_list]
    if args.lesser_than is not None:
        inlet_list = [take_lesser_than(inlet, args.lesser_than) for inlet in inlet_list]

    for inlet in inlet_list:
        print(inlet.name, "\n")
        print("Temperature Data:", inlet.temperature_data)
        print("Salinity Data:", inlet.salinity_data)
        print("Oxygen Data:", inlet.oxygen_data)
        print()


if __name__ == "__main__":
    main()