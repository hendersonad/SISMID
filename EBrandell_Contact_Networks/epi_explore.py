##### transmissibility vs epidemic size 

urb = open('urban_smaller.csv')
V = Graph()

for edge in urb:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    V.add_edge(node1,node2)
draw(V)