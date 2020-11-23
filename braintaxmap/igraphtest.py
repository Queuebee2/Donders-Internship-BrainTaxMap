from igraph import Graph, plot, layout;
import sys

g = Graph()
print('PREFIX:' ,sys.prefix)
print('INFO:', sys.version)

g.add_vertices(3)
g.add_edges([(0,1), (1,2)])

layout = g.layout("kk")
plot(g)
print(g)