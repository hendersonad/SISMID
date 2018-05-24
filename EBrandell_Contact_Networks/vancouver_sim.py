##### vancouver sim 

from networkx import *
from random import *
from pylab import hist, show

###### run the sim 100 times with vancouver ########
net = open('sexual_net.csv')
G = Graph()

for edge in net:
    node1,node2 = edge.strip().split(',') # split nodes by commas in CSV
    G.add_edge(node1,node2)
#draw(G)


T = 0.9  # probability of transmission

results = []
net_size = float(G.order())
  
for i in range(100):
    p_zero = choice(G.nodes())  # randomly choose a node
    infected = [p_zero]
    recovered = []
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
    print epi_size
    results.append(epi_size)

hist(results)
show()   # highly right skewed