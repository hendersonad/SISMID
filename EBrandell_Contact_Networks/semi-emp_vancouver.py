####### semi-empircal example ######

urb = open('urban_net.csv')
V = Graph()

for edge in urb:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    V.add_edge(node1,node2)
draw(V, with_labels=True)

dv = V.degree().values() # look at degree distribution
V.degree() # degree of each node
number_connected_components(V) # how many disjoint subgraphs = 1 (weird bc looks like 3 in graph)

V.order() # number of nodes, 2693
V.size() # number of edges, 22745
min(dv) # min degree = 2
max(dv) # max degree = 59
mean(dv) # mean degree = 16.8919

hist(dv) # one peak close to 15, ~N distribution