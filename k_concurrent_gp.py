import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math
import dwave
# from neal import SimulatedAnnealingSampler
import sys
#import dwave_networkx as dnx
import dimod


def solve_QUBO(G, nparts, vdegree):
    # ------- Set up our QUBO dictionary -------

    # Initialize our Q matrix
    Q = defaultdict(int)

    nNodes = len(G.nodes)
    penalty = math.ceil(min(2*max(vdegree),nNodes)/8)
    print("penalty ", penalty)
    alpha_lagr = 1000*[penalty]*nparts
    gamma_lagr = 1000*[penalty]*nNodes

    for i in range(nparts):
        for u in G.nodes:
            Q[(i*nNodes+u,i*nNodes+u)] += vdegree[u] + alpha_lagr[i] + gamma_lagr[u] -2*gamma_lagr[i] -2*(nNodes/nparts)*alpha_lagr[u]
        for v, w, d in G.edges(data=True):
            Q[(i*nNodes+v, i*nNodes+w)] += -d["weight"]
    # Run the QUBO on the solver from your config file

    sampler = EmbeddingComposite(DWaveSampler())

    # sampler = SimulatedAnnealingSampler()
    response = sampler.sample_qubo(Q,
                                   #chain_strength=chain_strength,
                                   num_reads=num_reads,
                                   label='Example - Graph Partitioning')

    sample = response.record.sample[0]

    return sample


def graphPart(G, nparts, vdegree):
    print_graph(G)

    sample = solve_QUBO(G, nparts, vdegree)

    return sample

def distance(G):
    distance=0
    for u, v in G.edges:
        distance += (abs(u-v)) * G[u][v]['weight']
    print("DISTANCE = {}.".format(distance))

def print_graph(G):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    print("nodes", G.nodes)
    print("edges", G.edges)

# ------- Set tunable parameters -------
num_reads = 1000
#gamma = 80

# ------- Set up our graph -------


def main():
    nq=8
    G = nx.Graph()
    G.add_nodes_from(range(nq))
    #G.add_weighted_edges_from([(0, 3, 1), (1, 2, 1), (1, 3, 2)])
    G.add_weighted_edges_from([(0, 3, 1), (1, 2, 1), (1, 3, 2), (2, 5, 3), (3, 7, 4), (4, 7, 1)])
    #G.add_weighted_edges_from([(3,4,1), (1,4,1),(0,5,1),(3,5,1),(0,6,1),(3,6,2),(2,5,1),(5,6,2),(2,6,1)]) # GRAPH FROM PAPER WITH SWAP COUNT 15
    
    distance(G)
    vdegree = list(dict(G.degree(weight="weight")).values())
    nparts = nq
    
    # RELABEL the graph to count the new distance
    total_ordering = graphPart(G, nparts, vdegree)
    mapping={}
    for i in range(len(total_ordering)):
        mapping[total_ordering[i]]=i

    print("total order")

    for i in range(nq):
        for j in range(nparts):
            print(total_ordering[i*nq+j], end=" ")
        print("\n")
    
    Gnew=nx.relabel_nodes(G,mapping)
    
    print_graph(Gnew)
    distance(Gnew)
    
    return 0


try:
    main()
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)
