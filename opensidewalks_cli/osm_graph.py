import networkx as nx
import osmium
import pyproj
from shapely.geometry import LineString


STREET_TAGS = [
    "primary",
    "secondary",
    "tertiary",
    "residential",
]


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

                # Only one exception: we shouldn't remove a node that's shared
                # between two different ways: this is an important decision
                # point for some paths.
                if edge_in["osm_id"] != edge_out["osm_id"]:
                    continue

                node_data = (node_in, node, node_out, edge_in["segment"])

                # Group by way
                edge_id = edge_in["osm_id"]
                if edge_id in remove_nodes:
                    remove_nodes[edge_id].append(node_data)
                else:
                    remove_nodes[edge_id] = [node_data]

        # NOTE: an otherwise unconnected circular path would be removed, as all
        # nodes are degree 2 and on the same way. This path is pointless for a
        # network, but is something to keep in mind for any downstream
        # analysis.
        for way_id, node_data in remove_nodes.items():
            # Sort by segment number
            sorted_node_data = list(sorted(node_data, key=lambda x: x[3]))

            # Group by neighboring segments
            groups = {}
            last_s = -10
            for ni, n, no, s in sorted_node_data:
                if (s - last_s) != 1:
                    # The last segment and this segment are not neighbors -
                    # create new group
                    receiving_edge = (ni, n)
                    groups[receiving_edge] = []
                groups[receiving_edge].append(n)
                last_s = s

            # Remove internal nodes by group
            # FIXME: ended up with two edges (!) in many cases. DUPES!
            for (u, v), nodes in groups.items():
                edge_data = self.G[u][v][0]
                ndref = edge_data["ndref"]
                for node in nodes:
                    # Append following to ndref
                    following_node = next(self.G.successors(node))
                    ndref.append(following_node)
                    self.G.remove_edge(node, following_node)
                self.G.add_edges_from([(u, node, edge_data)])

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
