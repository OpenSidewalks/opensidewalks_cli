import networkx as nx
import osmium
import pyproj
from shapely.geometry import LineString


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

        return OSMGraph(G)

    def simplify(self):
        """Simplifies graph by merging way segments of degree 2 - i.e.
        continuations.

        """
        # Structure is way_id: (node, segment_number). This makes it easy to
        # sort on-the-fly.
        remove_nodes = {}

        for node in self.G.nodes:
            if self._filter_nodes_to_remove(node):
                (
                    node_in,
                    node,
                    node_out,
                    incoming_way_id,
                    incoming_way_segment,
                ) = self._create_incoming_segment(node)

                # Group by way
                node_data = (node_in, node, node_out, incoming_way_segment)
                if incoming_way_id in remove_nodes:
                    remove_nodes[incoming_way_id].append(node_data)
                else:
                    remove_nodes[incoming_way_id] = [node_data]

        # NOTE: an otherwise unconnected circular path would be removed, as all
        # nodes are degree 2 and on the same way. This path is pointless for a
        # network, but is something to keep in mind for any downstream
        # analysis.
        for way_id, node_data in remove_nodes.items():
            self._merge_ways(node_data)

    def construct_geometries(self):
        # Given the current list of node references per edge, construct
        # geometry
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

    def _filter_nodes_to_remove(self, node):
        predecessors = list(self.G.predecessors(node))
        successors = list(self.G.successors(node))

        if (len(predecessors) == 1) and (len(successors) == 1):
            # Only one predecessor and one successor - ideal internal node
            # to remove from graph, merging its location data into other
            # edges.
            node_in = predecessors[0]
            node_out = successors[0]
            edge_in = self.G[node_in][node][0]
            edge_out = self.G[node][node_out][0]

            # Don't remove nodes shared by two different ways (for now).
            # FIXME: do we really want to do this?
            if edge_in["osm_id"] != edge_out["osm_id"]:
                return False
            return True
        return False

    def _create_incoming_segment(self, node):
        predecessors = list(self.G.predecessors(node))
        successors = list(self.G.successors(node))
        node_in = predecessors[0]
        node_out = successors[0]
        incoming_way = self.G[node_in][node][0]
        incoming_way_id = incoming_way["osm_id"]
        incoming_way_segment = incoming_way["segment"]

        return (node_in, node, node_out, incoming_way_id, incoming_way_segment)

    def _merge_ways(self, node_data):
        # node_data is a list of nodes to remove along with metadata about them
        # like the two neighboring nodes (ordered) and the index of the segment
        # preceding it.

        # Sort by segment index so that nodes are in same order as in the way
        sorted_node_data = list(sorted(node_data, key=lambda x: x[3]))

        # Group by neighboring segments
        groups = []

        node_in, node, node_out, segment_index = sorted_node_data[0]
        group = {
            "u": node_in,
            "v": node_out,
            "v_original": node,
            "merge_nodes": [node],
            "drop_edges": [(node, node_out)],
        }
        last_segment_index = segment_index
        for node_in, node, node_out, segment_index in sorted_node_data[1:]:
            if (segment_index - last_segment_index) == 1:
                # They're neighbors, so extend the current group
                group["v"] = node_out
                group["merge_nodes"].append(node)
                group["drop_edges"].append((node, node_out))
            else:
                # They aren't neighbors, so create a new group
                groups.append(group)
                group = {
                    "u": node_in,
                    "v": node_out,
                    "v_original": node,
                    "merge_nodes": [node],
                    "drop_edges": [(node, node_out)],
                }
            last_segment_index = segment_index

        # Append the last group
        groups.append(group)

        # Remove internal nodes by group
        for group in groups:
            u = group["u"]
            v = group["v"]
            v_original = group["v_original"]
            merge_nodes = group["merge_nodes"]

            edge_data = self.G[u][v_original][0]

            ndref = edge_data["ndref"]
            for merge_node in merge_nodes:
                # Add to node references
                ndref.append(merge_node)
            ndref.append(v)

            # Remove edges that have been merged.
            for u_drop, v_drop in group["drop_edges"]:
                self.G.remove_edge(u_drop, v_drop)

            self.G.add_edges_from([(u, v, edge_data)])
