# `opensidewalks_cli`

A cli application for generating OpenSidewalks task polygons for use with the
[HOT Tasking Manager](https://github.com/hotosm/tasking-manager). The task
polygons generated by this tool are optimized for two workflows:

- Pedestrian intersections, where a task will be centered on a street
intersection. This guarantees that mapping a crosswalk won't be interrupted by
crossing a task boundary.

- Sidewalks, where a task will encompass a single "block". If a "block" has
full sidewalk coverage, this task would be an outline around the entire
sidewalk path.

## Installation

### 1. Clone the repository

First, clone this repository: `git clone https://github.com/opensidewalks/opensidewalks_cli`

### 2. Installing the application

This tool is developed using [`poetry`](https://python-poetry.org/).

It can be installed with `poetry` by running `poetry install` in the main repo
directory or, if you have a recent version of `pip`, by running `pip install .`.

Alternative: Use the included Dockerfile to automatically install and build all
dependencies:

    cd <CLONED_DIR> && docker build -t opensidewalks_cli .

Then, for all subsequent commands that use `osw`, instead prepend with:
    docker run -v $(pwd)/path/to/data_directory:/data opensidewalks_cli


## Creating task polygons for intersection and sidewalk mapping

Run `osw task <path-to-osm-pbf> <path-to-aoi-polygon-geojson> <output-path>`,
where:

- `<path-to-osm-pbf` is a path to an OSM PBF file that has street data for the
area for which you want to create tasks.

- `<path-to-aoi-polygon-geojson>` is a path to a GeoJSON file with a single
polygon geometry defining your AOI (area of interest) - the boundary of your
tasks.

- `<output-path>` is your output directory, where several output files will be
created.

## Tutorial

We will use the example of Columbus, Ohio.

### 0. Prepare a work directory

Create and enter an empty directory:

`mkdir ohio_tasks && cd ohio_tasks`

Create subdirectories for our data:

`mkdir -p ./{data_sources,intermediate_data,output}`

### 1. Retrieve an OSM PBF

Download [OSM data for Ohio from geofabrik](http://download.geofabrik.de/north-america/us/ohio-latest.osm.pbf).

Create a directory called `data_sources` and place the OSM PBF in it.

`mv ~/Dowloads/ohio-latest.osm.pbf ./data_sources`

### 2. Define an area of interest polygon and store in a GeoJSON file

For this example, we will download and extract a polygon representing the North
Linden neighborhood of Columbus, Ohio. You can use any polygon that is of
interest to your project.

Download the [neighborhoods shapefile](https://opendata.columbus.gov/datasets/columbus-communities)
from Columbus Open Data Portal.

Place it in `data_sources`:

`mv ~/Downloads/Columbus_Commnities-shp.zip ./data_sources`. (Extract the zip
in this directory).

Open neighborhoods shapefile in QGIS, select 'North Linden' (objectid 383, area
name 'North Linden'), save a GeoJSON file with single polygon feature. Save to
the path `./intermediate_data/north_linden.geojson`.

(Optionally) Clip the OSM PBF using the AOI polygon and the osmium tool:
`osmium extract -p intermediate_data/north_linden.geojson data_sources/ohio-latest.osm.pbf -o intermediate_data/north_linden.osm.pbf`
. This will speed up the next steps if the extract is much larger than the
AOI.

### 3. Run the `osw` command line interface

Run the `task` subcommand of the `osw` cli application to create pedestrian
intersection and sidewalk task polygons.

`osw task intermediate_data/north_linden.osm.pbf intermediate_data/north_linden.geojson output/north_linden`

If you decided not to clip the OSM PBF, replace
`intermediate_data/north_linden.osm.pbf` with
`data_sources/ohio-latest.osm.pbf`.

If you receive an error that the `osw` command can't be found, make sure you
have activated the Python environment in which you installed
`opensidewalks_cli`. For example, if you cloned this repository and ran
`poetry install`, you will need to run `poetry shell` in this repository before
you can use `osw`. Alternatively, you can run the `osw` command from the base
directory of the repository by using `poetry run`, such as `poetry run osw`.
