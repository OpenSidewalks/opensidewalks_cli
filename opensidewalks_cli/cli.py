"""sidewalkscore CLI."""
import os

import click

from opensidewalks_cli.extract_tasks import (
    extract_crossing_tasks,
    extract_sidewalk_tasks,
)


@click.group()
def osw():
    pass


@osw.command()
@click.argument("osm_file", type=click.Path("r"))
@click.argument("area_geojson", type=click.Path("r"))
@click.argument("output_dir", type=click.Path("r"))
def task(osm_file, area_geojson, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    extract_crossing_tasks(osm_file, area_geojson, output_dir)
    extract_sidewalk_tasks(osm_file, area_geojson, output_dir)
