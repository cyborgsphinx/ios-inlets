Data Plotting Methods
=====================

The code in this repository has two main parts:

1. Aggregate and process data from multiple sources
1. Plot data in numerous configurations

The problem of sifting through the available data and producing plots was broken
down into many parts:

Data Aggregation
----------------

1. Read the data files
1. Check that a given file corresponds to a known inlet
1. Perform any necessary conversions or corrections
1. Categorize data based on depth

There are three main data file formats used: netCDF, IOS' internal file format,
and CSV. Each format is handled in a manner consistent with its structure, but
the subsequent steps remain the same for each.

Before fully processing the data, it is important to make sure that it is from
an inlet of interest. Each inlet is defined using a geoJSON polygon, allowing
for simple bounds checks to determine if it is of interest for any given inlet.

Some files contain data presented with different units than the ones displayed
in the plots, and must be converted. In these cases, a unit conversion method is
defined and used for all known unit conversion issues, using TEOS-10 functions
as needed.

As a final step in the aggregation process, the data is tagged based on its
depth. geoJSON allows for additional properties to be stored alongside the
geometry, and the particular boundaries for the depth categories are stored
there and used here for this purpose.

Plotting
--------

Plotting is done using [Matplotlib](https://matplotlib.org). Plots are broken
out into different kinds:

- Individual inlet data plots
- Individual inlet averaged data plots
- Yearly sampling histograms
- Monthly sampling heatmaps
- Annual averaging plots
- Annual anomaly plots

When called for, averaging is done in multiple stages in order to avoid bias
toward sources with more samples:

1. Data is averaged across a single day for a given file
1. These averages are combined to produce an average across a single month,
combining all data sources
1. If looking at yearly plots, the monthly averages are combined to produce a
yearly average
