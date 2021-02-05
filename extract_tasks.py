import json
import os
import sys

import geopandas as gpd
import networkx as nx
import numpy as np
import osmium
import pyproj
import scipy
from shapely.geometry import LineString, Polygon
from shapely.ops import polygonize


STREET_TAGS = [
    "primary",
    "secondary",
    "tertiary",
    "residential",
]


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
                return (
                    LineString([(point.x, point.y)] + coords[i + 1 :]),
                    coords[:i] + [(point.x, point.y)],
                )
            # print(pd, distance)

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
    lines.append(LineString(coords))
    return lines


def point_distance(p1, p2):
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    return (dx ** 2 + dy ** 2) ** 0.5


class OSMGraph:
    def __init__(self, G=None):
        if G is not None:
            self.G = G

        # Geodesic distance calculator. Assumes WGS84-like geometries.
        self.geod = pyproj.Geod(ellps="WGS84")

    @classmethod
    def from_pbf(self, pbf, way_filter=None):
        class OSMParser(osmium.SimpleHandler):
            def __init__(self):
                osmium.SimpleHandler.__init__(self)
                self.G = nx.MultiDiGraph()
                if way_filter is None:
                    self.way_filter = lambda u, v, d: True
                else:
                    self.way_filter = way_filter

            def way(self, w):
                if not self.way_filter(w.tags):
                    return

                d = {"osm_id": w.id, "tags": {}}

                for i in range(len(w.nodes) - 1):
                    u = w.nodes[i]
                    v = w.nodes[i + 1]
                    d2 = {**d}
                    d2["segment"] = i
                    d2["ndref"] = [u.ref, v.ref]
                    self.G.add_edges_from([(u.ref, v.ref, d2)])
                    self.G.add_node(u.ref, lon=u.lon, lat=u.lat)
                    self.G.add_node(v.ref, lon=v.lon, lat=v.lat)

        parser = OSMParser()
        parser.apply_file(pbf, locations=True)

        G = parser.G

        del parser

        # # TODO: allow selection of graph container type - in-memory vs.
        # entwiner
        # G = nx.MultiDiGraph()
        # # TODO: logging?

        # if way_filter is None:
        #     way_filter = lambda u, v, d: True

        # # Have to do two passes:
        # # Pass 1: extract graph structure using ways and node refs
        # # Pass 2: extract node lon-lat information so that we only retain
        # # used nodes
        # # in memory
        # # NOTE: if data structure was already ordered with ways first, then
        # # nodes, this
        # # would not require two passes. Way to verify?
        # osm_iter = osmread.parse_file(pbf)

        # print("Adding ways...")
        # # TODO: pre-filter incoming data with faster tool. Separate out nodes
        # # and ways,
        # # ideally with tag filter. Read as separate files. Lots of time
        # # wasted
        # # iterating over elements we don't analyze, because nodes occur
        # # before ways
        # # in these files and we need to know ways before nodes.
        # for primitive in osm_iter:
        #     if isinstance(primitive, osmread.Way):
        #         # Why does way_filter take u and v arguments anyways?
        #         if not way_filter(primitive.tags):
        #             continue

        #         d = {
        #             "osm_id": primitive.id,
        #             "tags": {}
        #         }
        #         d["ignore_incline"] = 0
        #         for key, values in IGNORE_INCLINE_TAGS.items():
        #             if key in primitive.tags and primitive.tags[key] in values:
        #                 d["ignore_incline"] = 1
        #         for key, values in KEEP_KEYS.items():
        #             if key in primitive.tags and primitive.tags[key] in values:
        #                 d["tags"][key] = primitive.tags[key]

        #         for i, (u, v) in enumerate(zip(primitive.nodes, primitive.nodes[1:])):
        #             d2 = {**d}
        #             d2["segment"] = i
        #             d2["ndref"] = [u, v]
        #             G.add_edges_from([(u, v, d2)])

        # print("Adding nodes...")
        # osm_iter = osmread.parse_file(pbf)
        # for primitive in osm_iter:
        #     if isinstance(primitive, osmread.Node):
        #         if primitive.id in G:
        #             # Add any relevant node data
        #             G.add_node(
        #                 primitive.id,
        #                 lon=round(primitive.lon, 7),
        #                 lat=round(primitive.lat, 7)
        #             )

        return OSMGraph(G)

    def simplify(self):
        """Simplifies graph by merging way segments of degree 2 - i.e. continuations."""
        # Structure is way_id: (node, segment_number). This makes it easy to sort
        # on-the-fly.
        remove_nodes = {}

        for node in self.G.nodes:
            predecessors = list(self.G.predecessors(node))
            successors = list(self.G.successors(node))

            if (len(predecessors) == 1) and (len(successors) == 1):
                # Only one predecessor and one successor - ideal internal node to
                # remove from graph, merging its location data into other edges.
                node_in = predecessors[0]
                node_out = successors[0]
                edge_in = self.G[node_in][node][0]
                edge_out = self.G[node][node_out][0]

                # Only one exception: we shouldn't remove a node that's shared between
                # two different ways: this is an important decision point for some
                # paths.
                if edge_in["osm_id"] != edge_out["osm_id"]:
                    continue

                node_data = (node_in, node, node_out, edge_in["segment"])

                # Group by way
                edge_id = edge_in["osm_id"]
                if edge_id in remove_nodes:
                    remove_nodes[edge_id].append(node_data)
                else:
                    remove_nodes[edge_id] = [node_data]

        # NOTE: an otherwise unconnected circular path would be removed, as all nodes
        # are degree 2 and on the same way. This path is pointless for a network, but
        # is something to keep in mind for any downstream analysis.
        for way_id, node_data in remove_nodes.items():
            # Sort by segment number
            sorted_node_data = list(sorted(node_data, key=lambda x: x[3]))

            # Group by neighboring segments
            groups = {}
            last_s = -10
            for ni, n, no, s in sorted_node_data:
                if (s - last_s) != 1:
                    # The last segment and this segment are not neighbors - create
                    # new group
                    receiving_edge = (ni, n)
                    groups[receiving_edge] = []
                groups[receiving_edge].append(n)
                last_s = s

            # Remove internal nodes by group
            # FIXME: ended up with two edges (!) in many cases. DUPES!
            for (u, v), nodes in groups.items():
                edge_data = self.G[u][v][0]
                ndref = edge_data["ndref"]
                last_node = u
                for node in nodes:
                    # Append following to ndref
                    following_node = next(self.G.successors(node))
                    ndref.append(following_node)
                    self.G.remove_edge(node, following_node)
                    last_node = node
                self.G.add_edges_from([(u, node, edge_data)])

    def construct_geometries(self):
        # Given the current list of node references per edge, construct geometry
        for u, v, d in self.G.edges(data=True):
            coords = []
            for ref in d["ndref"]:
                # FIXME: is this the best way to retrieve node attributes?
                node_d = self.G._node[ref]
                coords.append((node_d["lon"], node_d["lat"]))

            geometry = LineString(coords)
            d["_geometry"] = geometry
            d["length"] = self.geod.geometry_length(geometry)

        # FIXME: remove orphaned nodes!

    def to_undirected(self):
        if self.G.is_multigraph():
            G = nx.MultiGraph(self.G)
        else:
            G = nx.Graph(self.G)
        return OSMGraph(G)

    def get_graph(self):
        return self.G

    def filter_edges(self, func):
        # TODO: put this in a "copy-like" function
        if self.G.is_multigraph():
            if self.G.is_directed():
                G = nx.MultiDiGraph()
            else:
                G = nx.MultiGraph()
        else:
            if self.G.is_directed():
                G = nx.DiGraph()
            else:
                G = nx.Graph()

        for u, v, d in self.G.edges(data=True):
            if func(u, v, d):
                G.add_edge(u, v, **d)

        # Copy in node data
        for node in G.nodes:
            d = self.G._node[node]
            G.add_node(node, **d)

        return OSMGraph(G)

    def is_multigraph(self):
        return self.G.is_multigraph()

    def is_directed(self):
        return self.G.is_directed()


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
    gdf_area = gpd.GeoDataFrame(
        [{"geometry": Polygon(polygon["geometry"]["coordinates"][0])}]
    )
    gdf_trimmed = gpd.overlay(gdf_tasks, gdf_area, how="intersection")

    tasks_fc = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [list(row["geometry"].exterior.coords)],
                },
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


def extract_sidewalk_tasks(osm_file, neighborhood_geojson, output_dir):
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
    polygon_shapely = Polygon(polygon["geometry"]["coordinates"][0])
    area_line = LineString(polygon["geometry"]["coordinates"][0])
    gdf_polygon = gpd.GeoDataFrame(geometry=[polygon_shapely])
    # gdf_clip_points = gpd.overlay(gdf_lines, gdf_polygon, how="intersection")
    # print(gdf_clip_points)

    # print(
    #     LineString([[0, 0], [1, 1]]).intersection(LineString([[0, 1], [1, 0]]))
    # )

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

    polygon_lines = cut_polygon(area_line, distances, intersecting_points)

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
            "geometry": {
                "type": "Polygon",
                "coordinates": [list(polygon.exterior.coords)],
            },
            "properties": {},
        }
        features.append(feature)
    tasks_fc = {"type": "FeatureCollection", "features": features}

    output_path = os.path.join(output_dir, "sidewalk_tasks.geojson")
    with open(output_path, "w") as f:
        json.dump(tasks_fc, f)

    # print(next(iter(G.edges(data=True))))


if __name__ == "__main__":
    osm_file = sys.argv[1]
    neighborhood_geojson = sys.argv[2]
    output_dir = sys.argv[3]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    extract_crossing_tasks(osm_file, neighborhood_geojson, output_dir)
    extract_sidewalk_tasks(osm_file, neighborhood_geojson, output_dir)
