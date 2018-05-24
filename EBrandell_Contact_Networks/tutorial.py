file_in = open('simple_edges.txt')

for edge in file_in:
    print(edge)
    
file_in.close()

### first, make CSV from text file

file_in = open('simple_edges.txt')
file_out = open('simple_edges.csv', 'w')

for edge in file_in:
    node1,node2 = edge.split(' ')
    csv_edge = node1+','+node2
    file_out.write(csv_edge)
    
file_in.close()
file_out.close()

## make a network from scratch for example

from networkx import *
G = Graph()
G.add_edge('A','B')
G.add_edge('A','C')
G.add_node('E')
G.add_edge('E','C')
G.add_node('D')
G.add_edge('B','D')
G.add_edge('B','C')
G.add_node('F')
G.add_edge('C','F')

from pylab import show
draw(G)

G.remove_edge('A','C')
G.remove_node('A')
from pylab import show
draw(G, with_labels=True)

G.nodes()  # list of nodes
G.edges()  # list of edges
G.order()  # number of nodes
G.size()  # number of edges







