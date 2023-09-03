# Copyright 2019 D-Wave Systems, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ------ Import necessary packages ----
import networkx as nx
from collections import defaultdict
from itertools import combinations
#from dwave.system.samplers import DWaveSampler
#from dwave.system.composites import EmbeddingComposite
import math
import dwave
from neal import SimulatedAnnealingSampler
import sys
#import dwave_networkx as dnx
import dimod


def solve_QUBO(G):
    # ------- Set up our QUBO dictionary -------

    # Initialize our Q matrix
    Q = defaultdict(int)

    # Fill in Q matrix
    for u, v in G.edges:
        Q[(u, u)] += 1
        Q[(v, v)] += 1
        Q[(u, v)] += -2

    for i in G.nodes:
        Q[(i, i)] += gamma*(1-len(G.nodes))

    for i, j in combinations(G.nodes, 2):
        Q[(i, j)] += 2*gamma

    # ------- Run our QUBO on the QPU -------

    # Set chain strength
    chain_strength = gamma*len(G.nodes)

    # Run the QUBO on the solver from your config file

    # sampler = EmbeddingComposite(SimulatedAnnealingSampler())

    sampler = SimulatedAnnealingSampler()
    response = sampler.sample_qubo(Q,
                                   chain_strength=chain_strength,
                                   num_reads=num_reads,
                                   label='Example - Graph Partitioning')

    # See if the best solution found is feasible, and if so print the number of cut edges.
    sample = response.record.sample[0]
    # dwave.inspector.show(sample)
    print(sample)

    return sample


def graphPart(G, nparts, total_ordering):
    print("Graph on {} nodes created with {} out of {} possible edges.".format(
        len(G.nodes), len(G.edges), len(G.nodes) * (len(G.nodes)-1) / 2))
    print("nodes", G.nodes)
    print("edges", G.edges)

    print("nparts", nparts)
    if (nparts > 1):
        print("nparts > 1")

        sample = solve_QUBO(G)

        # In the case when n is odd, the set may have one more or one fewer nodes
        if sum(sample) in [math.floor(len(G.nodes)/2), math.ceil(len(G.nodes)/2)]:
            print("Valid partition found.")
        else:
            print("Invalid partition.")

        # partition graph into 2 parts
        lG = nx.Graph()
        rG = nx.Graph()
        for u in G.nodes:
            if sample[list(G.nodes).index(u)] == 0:
                lG.add_node(u)
            if sample[list(G.nodes).index(u)] == 1:
                rG.add_node(u)
        for u, v in G.edges:
            print("u, v", u, v)
            if sample[list(G.nodes).index(u)] == 0 and sample[list(G.nodes).index(v)] == 0:
                lG.add_edge(u, v)
            if sample[list(G.nodes).index(u)] == 1 and sample[list(G.nodes).index(v)] == 1:
                rG.add_edge(u, v)

        if (len(list(lG.nodes)) == 1 & len(list(rG.nodes)) == 1):
            print("2 nodes", list(lG.nodes) + list(rG.nodes))
            return list(lG.nodes) + list(rG.nodes)

        if (len(list(lG.nodes())) >= 2):
            total_ord_from_part_lG = graphPart(
                lG, math.floor(nparts/2), total_ordering)
            print("partial ordering AFTER PARTITIONING lG ",
                  total_ord_from_part_lG)
            print("total ordering before ", total_ordering)
            total_ordering = total_ordering + total_ord_from_part_lG
            print("total ordering after ", total_ordering)
        if (len(list(rG.nodes())) >= 2):
            total_ord_from_part_rG = graphPart(rG, math.floor(nparts/2), [])
            print("partial ordering AFTER PARTITIONING rG ",
                  total_ord_from_part_rG)
            print("total ordering before ", total_ordering)
            total_ordering = total_ordering + total_ord_from_part_rG
            print("total ordering after ", total_ordering)
    else:
        print("nparts <= 1")
    return total_ordering


# ------- Set tunable parameters -------
num_reads = 1000
gamma = 80

# ------- Set up our graph -------


def main():
    G = nx.Graph()
   # nq=8
    nq = 4
    for i in range(nq):
        G.add_node(i)
   # G.add_edges_from([(0, 5), (1, 6), (1, 7), (2, 5), (3, 7), (4, 7)])
    G.add_edges_from([(0,3), (1,3), (0,2)])
    nparts = nq
    total_ordering = graphPart(G, nparts, [])
    print("total order", total_ordering)

    return 0


try:
    main()
except:
    import traceback
    traceback.print_exc()
    sys.exit(1)
