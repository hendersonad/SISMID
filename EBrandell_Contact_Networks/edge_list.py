###### AUTOMATING NETWORK CONSTRUCTION #####

edgelist = open('toy_edges.csv')
for edge in edgelist:
    print edge
edgelist.close()

edgelist = open('toy_edges.csv')

G = Graph()

for edge in edgelist:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    G.add_edge(node1,node2)
draw(G, with_labels=True)

edgelist.close()
