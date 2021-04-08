import json
import os

import geopandas as gpd
import numpy as np
import scipy
from shapely.geometry import (
    MultiPolygon,
    Polygon,
    mapping,
    shape,
)

from opensidewalks_cli.exceptions import InvalidPolygonError
from opensidewalks_cli.osm_graph import OSMGraph


STREET_TAGS = [
    "primary",
    "secondary",
    "tertiary",
    "residential",
]


def extract_crossing_tasks(osm_file, neighborhood_geojson, output_dir):
    # Load neighborhood
    with open(neighborhood_geojson) as f:
        polygon_fc = json.load(f)
    polygon = polygon_fc["features"][0]

    # Loop through the OSM file and extract network + coordinate information
    def way_filter(d):
        if "highway" not in d:
            return False
        if d["highway"] in STREET_TAGS:
            return True
        return False

    OG = OSMGraph.from_pbf(osm_file, way_filter=way_filter)
    OG.simplify()
    G = OG.G
    # Extract intersections
    nodes = list(G.nodes)
    for node in nodes:
        if G.degree(node) <= 2:
            del G._node[node]

    points = [(d["lon"], d["lat"]) for n, d in G.nodes(data=True)]
    # polygon_coords = polygon["geometry"]["coordinates"][0][:-1]
    # points += polygon_coords
    points = np.array(points)

    points_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": list(point)},
                "properties": {},
            }
            for point in points
        ],
    }
    output_path = os.path.join(output_dir, "crossing_tasks_points.geojson")
    with open(output_path, "w") as f:
        json.dump(points_fc, f)

    # Create voronoi polygons
    vor = scipy.spatial.Voronoi(points)
    features = []
    for region in vor.regions:
        if not region:
            continue
        if -1 in region:
            continue
        coords = [vor.vertices[index] for index in region]
        # Return to first coordinate to create valid polygon
        coords.append(coords[0])
        feature = {"geometry": Polygon(coords)}
        features.append(feature)

    # Trim task polygons
    gdf_tasks = gpd.GeoDataFrame(features)

    try:
        polygon_shapely = shape(polygon["geometry"])
    except ValueError:
        raise InvalidPolygonError()

    if polygon_shapely.type not in ["Polygon", "MultiPolygon"]:
        raise InvalidPolygonError()

    if polygon_shapely.type == "Polygon":
        polygon_shapely = MultiPolygon([polygon_shapely])

    gdf_area = gpd.GeoDataFrame([{"geometry": polygon_shapely}])
    gdf_trimmed = gpd.overlay(gdf_tasks, gdf_area, how="intersection")

    features = []
    tasks_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(row["geometry"]),
                "properties": {},
            }
            for idx, row in gdf_trimmed.iterrows()
        ],
    }

    # FIXME: there are point without voronoi polygons near the edges.
    # Troubleshoot.

    # Write to disk
    output_path = os.path.join(output_dir, "crossing_tasks.geojson")
    with open(output_path, "w") as f:
        json.dump(tasks_fc, f)
