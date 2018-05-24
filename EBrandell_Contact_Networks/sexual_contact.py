##### sexual contact network ######
snet = open('sexual_net.csv')
S = Graph()

for edge in snet:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    S.add_edge(node1,node2)
draw(S, with_labels=True)

ds = S.degree().values() # look at degree distribution
S.degree() # degree of each node
number_connected_components(S) # how many disjoint subgraphs = 4

S.order() # number of nodes, 285
S.size() # number of edges, 287
min(ds) # min degree = 1
max(ds) # max degree = 9
mean(ds) # mean degree = 2.014

hist(ds) # one peak close to 1, right skewed (maybe neg binomial), 

