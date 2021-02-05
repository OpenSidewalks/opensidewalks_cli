"""sidewalkscore CLI."""
import os

import click

from opensidewalks_cli.exceptions import InvalidPolygonError
from opensidewalks_cli.extract_tasks import (
    extract_crossing_tasks,
    extract_sidewalk_tasks,
)


@click.group()
def osw():
    pass


@osw.command()
@click.argument(
    "osm_file", type=click.Path(file_okay=True, dir_okay=False, readable=True)
)
@click.argument(
    "area_geojson",
    type=click.Path(file_okay=True, dir_okay=False, readable=True),
)
@click.argument("output_dir", type=click.Path(exists=False))
def task(osm_file, area_geojson, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    try:
        extract_crossing_tasks(osm_file, area_geojson, output_dir)
        extract_sidewalk_tasks(osm_file, area_geojson, output_dir)
    except InvalidPolygonError:
        raise InvalidPolygonError(
            "Input polygon file is invalid - not a polygon."
        )
