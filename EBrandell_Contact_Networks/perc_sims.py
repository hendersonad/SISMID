####### percolation simulation #########

from networkx import *
from random import *

# define parameters
T = 0.9  # probability of transmission

snet = open('sexual_net.csv') # using sexual network to make network

# create network
G = Graph()

for edge in snet:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    G.add_edge(node1,node2)
#draw(G, with_labels=True)

net_size = float(G.order()) # 6

# node status
p_zero = choice(G.nodes())  # randomly choose a node
infected = [p_zero]
recovered = []

# run the sim
while len(infected) > 0:
    v = infected.pop(0)
    for u in G.neighbors(v):
        # keep node u susceptible
        if u not in infected and u not in recovered:
            if random() < T:
                infected.append(u)
    recovered.append(v)
    
# tally total fraction that gon infected
epi_size = len(recovered)/net_size
print epi_size  # yay! worked!












