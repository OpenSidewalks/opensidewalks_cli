import json
import os

import geopandas as gpd
from shapely.geometry import (
    LineString,
    MultiLineString,
    mapping,
    shape,
)
from shapely.ops import polygonize

from opensidewalks_cli.osm_graph import OSMGraph
from opensidewalks_cli.way_filters import street_filter


def cut_polygon(polygon, distances, points):
    # Move in the order of the distance along the polygon
    # Have to sort distances and points together
    zipped = zip(distances, points)
    zipped = sorted(zipped, key=lambda x: x[0], reverse=True)
    distances, points = zip(*zipped)

    coords = list(polygon.coords)

    def cut_line(distance, point, coords):
        pd = 0
        last = coords[0]
        for i, p in enumerate(coords[1:]):
            pd += point_distance(last, p)

            if pd > distance:
                j = i + 1
                return (
                    LineString([(point.x, point.y)] + coords[j:]),
                    coords[:i] + [(point.x, point.y)],
                )

            last = p
        # Should not be able to reach this situation - catch it
        return None, None
        print(pd, distance)
        raise Exception()
        # raise Exception()

    # Note: also cuts at the start/end point by default
    lines = []
    counter = 0
    for distance, point in zip(distances, points):
        new_line, new_coords = cut_line(distance, point, coords)
        if new_line is None:
            continue
        lines.append(new_line)
        coords = new_coords
        counter += 1
        if counter > 1:
            pass
            # break
        # pd = 0
        # last = coords[0]
        # for i, p in enumerate(coords[1:]):
        #     pd += point_distance(last, p)

        #     if pd > distance:
        #         lines.append(LineString([(point.x, point.y)] + coords[i:]))
        #         coords = coords[:i] + [(point.x, point.y)]
        #         break

        #     last = p
        # # This should only be reached in the case of some kind of error -
        # # catch it!
        # raise Exception()
        # # lines.append(LineString([(point.x, point.y)] + coords[i:]))
        # # coords = coords[:i] + [(point.x, point.y)]
    try:
        line = LineString(coords)
        lines.append(line)
    except ValueError:
        # TODO: implement logging, log as info
        pass
    return lines


def point_distance(p1, p2):
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    return (dx ** 2 + dy ** 2) ** 0.5


def extract_sidewalk_tasks(osm_file, neighborhood_geojson, output_dir):
    # Load neighborhood
    with open(neighborhood_geojson) as f:
        polygon_fc = json.load(f)
    polygon = polygon_fc["features"][0]

    OG = OSMGraph.from_pbf(osm_file, way_filter=street_filter)

    # OG.simplify()
    G = OG.G

    # Create list of lines from the graph
    lines = []
    line_features = []
    for u, v, d in G.edges(data=True):
        if d["ndref"]:
            coords = []
            for n in d["ndref"]:
                node_data = G.nodes[n]
                coords.append([node_data["lon"], node_data["lat"]])
            line = LineString(coords)
            lines.append(line)
            line_feature = {
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": coords},
                "properties": {"u": u, "v": v},
            }
            line_features.append(line_feature)

    line_fc = {
        "type": "FeatureCollection",
        "features": line_features,
    }
    output_path = os.path.join(output_dir, "sidewalk_tasks_lines.geojson")
    with open(output_path, "w") as f:
        json.dump(line_fc, f)

    # Need to deal with boundary issue - need to close on the polygon or
    # using a buffer. Two steps: (1) Trim back street lines based on
    # neighborhood polygon and (2) split neighborhood polygon into line
    # segments split up by where it meets the streets.

    # Create split-up neighborhoods polygon
    gdf_lines = gpd.GeoDataFrame(geometry=lines)
    polygon_shapely = shape(polygon["geometry"])
    if polygon_shapely.type == "Polygon":
        area_line = LineString(polygon_shapely.exterior.coords)
    else:
        area_line = MultiLineString(
            [geom.exterior.coords for geom in polygon_shapely]
        )
    gdf_polygon = gpd.GeoDataFrame(geometry=[polygon_shapely])

    gdf_intersecting_lines = gpd.sjoin(
        gdf_lines, gdf_polygon, op="intersects", how="inner"
    )
    distances = []
    intersecting_points = []
    for idx, geom in gdf_intersecting_lines["geometry"].iteritems():
        intersecting_geom = area_line.intersection(geom)
        if intersecting_geom.type == "LineString":
            continue
        if intersecting_geom.type == "MultiPoint":
            intersecting_geoms = list(intersecting_geom.geoms)
        else:
            intersecting_geoms = [intersecting_geom]

        for geom in intersecting_geoms:
            distances.append(area_line.project(geom))
            intersecting_points.append(geom)

    ixn_points_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [point.x, point.y],
                },
                "properties": {},
            }
            for point in intersecting_points
        ],
    }
    ixn_pt_path = os.path.join(output_dir, "sidewalk_tasks_ixn_points.geojson")
    with open(ixn_pt_path, "w") as f:
        json.dump(ixn_points_fc, f)

    if area_line.type == "LineString":
        polygon_lines = cut_polygon(area_line, distances, intersecting_points)
    else:
        for line in area_line:
            polygon_lines = []
            polygon_lines += cut_polygon(line, distances, intersecting_points)

    polygon_path = os.path.join(
        output_dir, "sidewalk_tasks_cut_polygon.geojson"
    )
    with open(polygon_path, "w") as f:
        polygon_fc = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": list(line.coords),
                    },
                    "properties": {},
                }
                for line in polygon_lines
            ],
        }
        json.dump(polygon_fc, f)

    # Trim the street lines
    gdf_trimmed = gpd.overlay(gdf_lines, gdf_polygon, how="intersection")

    # Combine the line data sources
    combined_lines = polygon_lines + list(gdf_trimmed["geometry"])

    # Write to disk
    sidewalk_tasks_polygons = polygonize(combined_lines)
    features = []
    for polygon in sidewalk_tasks_polygons:
        feature = {
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {},
        }
        features.append(feature)
    tasks_fc = {"type": "FeatureCollection", "features": features}

    output_path = os.path.join(output_dir, "sidewalk_tasks.geojson")
    with open(output_path, "w") as f:
        json.dump(tasks_fc, f)
