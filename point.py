# How to plot a point
import osmnx as ox
import matplotlib.pyplot as plt

G = ox.graph_from_address('2801 Debarr Road, Anchorage, AK 99508', distance=1600)

fig, ax = ox.plot_graph(G, show=False, close=False)

ax.scatter(-149.824945, 61.210608, c='red')
plt.show()

# pre 1971
# 1971-1990
# post 1990
