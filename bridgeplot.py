import osmnx as ox
import matplotlib.pyplot as plt
import xlrd

loc = ("bridges.xlsx")
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)

# get a graph for some city
G = ox.graph_from_address('3211 Providence Dr, Anchorage, AK 99508', distance=10000)

fig, ax = ox.plot_graph(G, show=False, close=False)

for i in range (1,130):
    y = sheet.row_values(i, start_colx=1, end_colx=2)
    x = sheet.col_values(2, start_rowx=i, end_rowx=i+1)
    [z] = sheet.col_values(3, start_rowx=i, end_rowx=i+1)
    if z <= 1971:
        ax.scatter(x,y, c='red')
    elif z > 1971 and z <= 1990:
        ax.scatter(x, y, c='yellow')
    else:
        ax.scatter(x,y, c='green')

plt.show()





