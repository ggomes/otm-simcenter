# Given three points off of the network, matches to the nearest edge and highlights that edge

import osmnx as ox
import matplotlib.pyplot as plt
ox.config(use_cache=True, log_console=True)

G = ox.graph_from_address('2801 Debarr Road, Anchorage, AK 99508', distance=400)

fig, ax = ox.plot_graph(G, show=False, close=False)

ax.scatter(-149.824945, 61.210608, c='red')
ax.scatter(-149.828, 61.2078, c='blue')
ax.scatter(-149.825, 61.2114, c='violet')
plt.show()

x = (-149.824945, -149.828, -149.825)
y = (61.210608, 61.2078, 61.2114)

edge = ox.get_nearest_edges(G, x, y, method='kdtree', dist=0.0001)

for i in range(0, len(edge)):
    ec = ['orange' if (u == edge[i, 1] and v == edge[i, 0]) else 'grey' for u, v, k in G.edges(keys=True)]
    fig, ax = ox.plot_graph(G, node_color='w', node_edgecolor='k', node_size=30,
                            node_zorder=3, edge_color=ec, edge_linewidth=2)



