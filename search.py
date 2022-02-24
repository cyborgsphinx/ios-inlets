import argparse
import inlets


def take_before(inlet, before):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.time.year < before
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.time.year < before
    ]
    inlet.oxygen_data = [
        datum for datum in inlet.oxygen_data if datum.time.year < before
    ]
    return inlet


def take_after(inlet, after):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.time.year > after
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.time.year > after
    ]
    inlet.oxygen_data = [
        datum for datum in inlet.oxygen_data if datum.time.year > after
    ]
    return inlet


def take_greater_than(inlet, greater_than):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.datum > greater_than
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.datum > greater_than
    ]
    inlet.oxygen_data = [
        datum for datum in inlet.oxygen_data if datum.datum > greater_than
    ]
    return inlet


def take_lesser_than(inlet, lesser_than):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.datum < lesser_than
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.datum < lesser_than
    ]
    inlet.oxygen_data = [
        datum for datum in inlet.oxygen_data if datum.datum < lesser_than
    ]
    return inlet


def take_deeper(inlet, depth):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.depth >= depth
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.depth >= depth
    ]
    inlet.oxygen_data = [datum for datum in inlet.oxygen_data if datum.depth >= depth]
    return inlet


def take_used(inlet):
    inlet.temperature_data = [
        datum for datum in inlet.temperature_data if datum.bucket != inlets.IGNORE
    ]
    inlet.salinity_data = [
        datum for datum in inlet.salinity_data if datum.bucket != inlets.IGNORE
    ]
    inlet.oxygen_data = [
        datum for datum in inlet.oxygen_data if datum.bucket != inlets.IGNORE
    ]
    return inlet


def print_data(inlet, **kwargs):
    only_filename = kwargs["filenameonly"] if "filenameonly" in kwargs else False
    use_data = kwargs["use_data"] if "use_data" in kwargs else "all"

    print(inlet.name, "\n")
    if only_filename:
        if use_data.lower() in ["all", "temperature"]:
            print(
                "Temperature Data:",
                [datum.filename for datum in inlet.temperature_data],
            )
        if use_data.lower() in ["all", "salinity"]:
            print("Salinity Data:", [datum.filename for datum in inlet.salinity_data])
        if use_data.lower() in ["all", "oxygen"]:
            print("Oxygen Data:", [datum.filename for datum in inlet.oxygen_data])
    else:
        if use_data.lower() in ["all", "temperature"]:
            print("Temperature Data:", inlet.temperature_data)
        if use_data.lower() in ["all", "salinity"]:
            print("Salinity Data:", inlet.salinity_data)
        if use_data.lower() in ["all", "oxygen"]:
            print("Oxygen Data:", inlet.oxygen_data)
    print()


def main():
    parser = argparse.ArgumentParser()
    # inlet retrieval args
    parser.add_argument("-r", "--from-saved", action="store_true")
    parser.add_argument("-n", "--skip-netcdf", action="store_true")
    parser.add_argument("-d", "--data", type=str, nargs="?", default="data")
    # search args
    parser.add_argument(
        "-b", "--before", metavar="YEAR", type=int, nargs="?", default=None
    )
    parser.add_argument(
        "-a", "--after", metavar="YEAR", type=int, nargs="?", default=None
    )
    parser.add_argument("-g", "--greater-than", type=float, nargs="?", default=None)
    parser.add_argument("-l", "--lesser-than", type=float, nargs="?", default=None)
    parser.add_argument("-i", "--inlet-name", type=str, nargs="?", default=None)
    parser.add_argument("-t", "--only-temperature", action="store_true")
    parser.add_argument("-s", "--only-salinity", action="store_true")
    parser.add_argument("-o", "--only-oxygen", action="store_true")
    parser.add_argument("-u", "--only-used", action="store_true")
    parser.add_argument("-m", "--depth", type=float, nargs="?", default=None)
    parser.add_argument("--file-name", action="store_true")
    args = parser.parse_args()
    inlet_list = inlets.get_inlets(args.data, args.from_saved, args.skip_netcdf)

    if args.before is not None:
        inlet_list = [take_before(inlet, args.before) for inlet in inlet_list]
    if args.after is not None:
        inlet_list = [take_after(inlet, args.after) for inlet in inlet_list]
    if args.greater_than is not None:
        inlet_list = [
            take_greater_than(inlet, args.greater_than) for inlet in inlet_list
        ]
    if args.lesser_than is not None:
        inlet_list = [take_lesser_than(inlet, args.lesser_than) for inlet in inlet_list]
    if args.depth is not None:
        inlet_list = [take_deeper(inlet, args.depth) for inlet in inlet_list]

    if any(
        [
            args.only_temperature and args.only_salinity,
            args.only_temperature and args.only_oxygen,
            args.only_salinity and args.only_oxygen,
        ]
    ):
        print(
            "--only-temperature, --only-salinity, and --only-oxygen are mutually exclusive"
        )
        exit(1)
    elif args.only_temperature:
        data = "temperature"
    elif args.only_salinity:
        data = "salinity"
    elif args.only_oxygen:
        data = "oxygen"
    else:
        data = "all"

    if args.only_used:
        inlet_list = [take_used(inlet) for inlet in inlet_list]

    if args.inlet_name is not None:
        inlet = [
            inlet
            for inlet in inlet_list
            if inlet.name.lower() == args.inlet_name.lower()
        ][0]
        print_data(inlet, filenameonly=args.file_name, use_data=data)
    else:
        for inlet in inlet_list:
            print_data(inlet, filenameonly=args.file_name, use_data=data)


if __name__ == "__main__":
    main()
