import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import osmnx as ox
from descartes import PolygonPatch
from shapely.geometry import Point, LineString, Polygon
ox.config(log_console=True, use_cache=True)
ox.__version__

# configure the place, network type, trip times, and travel speed
place = 'Anchorage, Alaska, USA'
network_type = 'drive'
trip_times = [10, 20, 30, 40, 50] # minutes
travel_speed = 50 # driving speed km/hour

# download street network
north, south, east, west = 61.2366, 61.1185, -149.6916, -149.9483
G = ox.graph_from_bbox(north, south, east, west, network_type='drive')

# find center most node and then project the graph to UTM
gdf_nodes = ox.graph_to_gdfs(G, edges=False)
x, y = gdf_nodes['geometry'].unary_union.centroid.xy
center_node = ox.get_nearest_node(G, (y[0], x[0]))
G = ox.project_graph(G)

# add an edge attribute for time in minutes required to traverse each edge
meters_per_minute = travel_speed * 100 / 60 # km/hr to m/min
for u, v, k, data in G.edges(data=True, keys=True):
    data['time'] = data['length'] / meters_per_minute


# get one color for each isochrone
iso_colors = ox.get_colors(n=len(trip_times), cmap='inferno', start=0.3, return_hex=True)

# color the nodes according to isochrone then plot the street network
node_colors = {}
for trip_time, color in zip(sorted(trip_times, reverse=True), iso_colors):
    subgraph = nx.ego_graph(G, center_node, radius=trip_time, distance='time')
    for node in subgraph.nodes():
            node_colors[node] = color
nc = [node_colors[node] if node in node_colors else 'none' for node in G.nodes()]
ns = [20 if node in node_colors else 0 for node in G.nodes()]
fig, ax = ox.plot_graph(G, fig_height=8, node_color=nc, node_size=ns, node_alpha=0.8, node_zorder=2)
