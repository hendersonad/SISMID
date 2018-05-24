#!/usr/bin/python
from networkx import *
from random import *
from pylab import mean

## note: running object. "tab" will give you list of commands you can execute on that object!

### Define percolation simulation, using graph G and transmissibility T
def percolate(G, T):
    states = dict([(node, 's') for node in G.nodes()]) 

    p_zero = choice(G.nodes())  # randomly infect an initial node
    states[p_zero] = 'i'
    infected = [p_zero]
    recovered = []  # empty retainer for recovered

    while len(infected) > 0:
        v = infected.pop(0)
        for u in G.neighbors(v):
            if states[u] == 's' and random() < T:  # infect if random number
                states[u] = 'i'
                infected.append(u)
        states[v] = 'r'
        recovered.append(v)
    ### return the epi size as a fraction of the ORIGINAL population
    return len(recovered)/float(net_size)


### Set transmissibility
### Corresponds to an R0 of about 2.5 for the urban network
T = 0.1357 

### Build urban network
file = open("urban_net.csv")
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()


############################################
### Implement intervention strategy here ###
############################################
###### run without vacc for comparison ######
file = open("urban_net.csv")
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()
results_s = []
for i in range(100):    ### Run the simulation a bunch of times
    s =  percolate(G, T)
    results_s.append(s)
### Print the mean epidemic size
print mean(results_s)
hist(results_s)

##### vaccinate 17.8% of population randomly by removing individuals with 0.178 chance #####
file = open("urban_net.csv")
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()

V = 0.178
results_s1 = []

for node in G.nodes():
    if random() < V:  # infect if random number
        G.remove_edges_from(G.edges(node)) # don't want to remove the node b/c that changes pop size             
for i in range(100):    ### Run the simulation a bunch of times
    s1 =  percolate(G, T)
    results_s1.append(s1)
### Print the mean epidemic size
print mean(results_s1)

hist(results_s1)



####### run vaccine that vaccinates top 17.8% degrees of populaiton #####
### Build urban network
file = open("urban_net.csv")
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()

results_s2 = []
for node in G.nodes():
    if G.degree(node) > 22:  # if degree is >=23, then it is top 17.8% of degree distribution
        G.remove_edges_from(G.edges(node)) # don't want to remove the node b/c that changes pop size             
for i in range(100):    ### Run the simulation a bunch of times
    s2 =  percolate(G, T)
    results_s2.append(s2)
### Print the mean epidemic size
print mean(results_s2)

hist(results_s2)


###### randomly reduce each node's degree by 17.8% #####
file = open("urban_net.csv")
G = Graph()
for edge in file:
    G.add_edge(*edge.strip().split(','))
net_size = G.order()

V = 0.178
results_s3 = []

for e in G.edges():
    if random() < V:  # infect if random number
        G.remove_edge(e[0],e[1]) # removal of edges from BOTH nodes involved 
for i in range(100):    ### Run the simulation a bunch of times
    s3 =  percolate(G, T)
    results_s3.append(s3)
### Print the mean epidemic size
print mean(results_s3)

hist(results_s3)

