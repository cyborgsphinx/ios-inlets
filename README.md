# Inlet Charts

This project is intended to update several inlet charts, as well as automate the process for the future.
The charts are visible at

- [Annual averaged deep water properties for inlets](https://www.pac.dfo-mpo.gc.ca/science/oceans/bc-inlets-mer-de-bras-cb/water-prop-eau-eng.html)
- [Long term trends in deep water properties of BC inlets](https://www.pac.dfo-mpo.gc.ca/science/oceans/bc-inlets-mer-de-bras-cb/index-eng.html)

The backing data can be found at [the waterproperties archive](https://www.waterproperties.ca/osd_data_archive/netCDF_Data/) and can be downloaded to `data/` for offline access:

    $ wget -m -np --cut-dirs=2 -P data https://www.waterproperties.ca/osd_data_archive/

Inlet polygons defined using https://geojson.io

## Dependencies

This project depends on (at least)

- python version 3.9 or later
- [poetry](https://python-poetry.org)
- a C compiler
- python development headers (python-dev or equivalent)

Once those are installed, you can run

    $ poetry install --no-dev

to get all the rest of the dependencies. To get dependencies which are useful for development, simply run

    $ poetry install

instead.

## Usage

To get temperature, salinity, and dissovled oxygen plots of all the inlets defined in `inlets.geojson`, run

    $ poetry run plot

To get monthly averages for temperature, salinity, and dissolved oxygen for all inlets, run

    $ poetry run plot -a

To get annual plots for all inlets, run

    $ poetry run plot -A

## GeoJSON Properties

When adding water bodies, certain property keys are picked up and added to the python object to influence its behaviour:

- "name": Used as an identifier for the data, in plot titles, file names, sqlite tables, and for filtering. Required.
    Example: `"name": "Saanich Inlet"`
- "area": Used to group inlets together for aggregate plots like the annual averages and annual anomalies charts. Required.
    Example: `"area": "Salish Sea"`
- "boundaries": Used to define the three depth categories for the deep water plots. Required.
    Example: `"boundaries": "[250, 350, 450]"`
- "shallow boundaries": Used to define the two depth categories for the shallow water plots, defaults to `[0, 30, 100]`.
    Example: `"shallow boundaries": "[0, 30, 100]"`
- "limits": Used to define view limits for T, S, O charts, broken out by "deep" and "shallow". An empty list disables the feature.
    Example: `"limits": { "oxygen": { "deep": [0, 5], "shallow": [] } }`
- "seasons": Used to define the seasons for seasonal trend removal charts. Expected to be the numbers `1..12` separated into a number of lists.
    Example: `"seasons": [[1, 2, 3], [4, 5], [6, 7, 8, 9], [10, 11, 12]]`
