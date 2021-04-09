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
from opensidewalks_cli.way_filters import street_filter

# TODO: cluster boulevards like OSMnx (or just use OSMnx for compatibility).

# TODO: should include intersections between streets and parking lot ways
#       but not intersections solely between parking lot ways


def extract_crossing_tasks(osm_file, neighborhood_geojson, output_dir):
    # Load neighborhood
    with open(neighborhood_geojson) as f:
        polygon_fc = json.load(f)
    polygon = polygon_fc["features"][0]

    OG = OSMGraph.from_pbf(osm_file, way_filter=street_filter)
    OG.simplify()
    G = OG.G

    # Save simplified graph data for debugging
    OG.construct_geometries()
    G_simplified_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(d["_geometry"]),
                "properties": {"u": u, "v": v},
            }
            for u, v, d in G.edges(data=True)
        ],
    }
    output_path = os.path.join(
        output_dir, "crossing_tasks_simplified_ways.geojson"
    )
    with open(output_path, "w") as f:
        json.dump(G_simplified_fc, f)

    # Extract intersections
    nodes = list(G.nodes)
    for node in nodes:
        if G.degree(node) <= 2:
            del G._node[node]

    point_data = [(n, (d["lon"], d["lat"])) for n, d in G.nodes(data=True)]
    # polygon_coords = polygon["geometry"]["coordinates"][0][:-1]
    # points += polygon_coords
    point_data = np.array(point_data)

    points_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": list(point)},
                "properties": {"node_id": n},
            }
            for n, point in point_data
        ],
    }
    output_path = os.path.join(output_dir, "crossing_tasks_points.geojson")
    with open(output_path, "w") as f:
        json.dump(points_fc, f)

    # Create voronoi polygons
    xy = np.array([p for n, p in point_data])
    if not xy.shape[0]:
        raise Exception(
            "No intersections! Is your area of interest within the OSM file?"
        )
    vor = scipy.spatial.Voronoi(xy)
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
