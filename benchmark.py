from src.serial.bellman_ford import bellman_ford
from src.serial.dijkstra import dijkstra
from src.serial.delta_stepping import delta_stepping
from src.serial.radius_stepping import radius_stepping

from src.parallel.delta_stepping_pl import delta_stepping_multiprocessing

import networkx as nx
import time
import random


NUM_VERTICES = 300


g = nx.complete_graph(NUM_VERTICES)
G = nx.Graph()
for u, v in g.edges:
    G.add_edge(u, v, weight=random.uniform(1, 10))

# G = nx.gnm_random_graph(15, 50)
# # G = nx.connected_watts_strogatz_graph(1000, 20, 0.5, tries=100, seed=1)
# for (u,v,data) in G.edges(data=True):
#     data['weight'] = round(random.uniform(1,10), 1)

if __name__ == '__main__':
    source = 0
    delta = 5
    r = 7

    print("SERIAL:\n")

    bf_start = time.perf_counter()
    bellman_ford_result = bellman_ford(G, source)
    bf_duration = time.perf_counter() - bf_start
    print(f"{'Bellman-Ford Time:':<25} {bf_duration:.20f}")

    dij_start = time.perf_counter()
    dijkstra_result = dijkstra(G, source)
    dij_duration = time.perf_counter() - dij_start
    print(f"{'Dijkstra Time:':<25} {dij_duration:.20f}")

    delta_step_start = time.perf_counter()
    delta_stepping_result = delta_stepping(G, source, delta)
    delta_duration = time.perf_counter() - delta_step_start
    print(f"{'Delta-stepping Time:':<25} {delta_duration:.20f}")

    radius_step_start = time.perf_counter()
    radius_stepping_result = radius_stepping(G, source, delta)
    radius_duration = time.perf_counter() - radius_step_start
    print(f"{'Radius-stepping Time:':<25} {radius_duration:.20f}")

    assert(delta_stepping_result == bellman_ford_result)
    assert(radius_stepping_result == bellman_ford_result)

    print("Correctness checks passed!\n")
    # print("Bellman-Ford Distances:", bellman_ford_result)
    # print("Dijkstra Distances:", dijkstra_result)
    # print("Delta-stepping Distances:", delta_stepping_result)
    # print("Radius-stepping Distances:", radius_stepping_result)

    print("PARALLEL:\n")

    delta_step_pl_start = time.perf_counter()
    delta_stepping_pl_result = delta_stepping_multiprocessing(G, source, delta)
    delta_pl_duration = time.perf_counter() - delta_step_pl_start
    print(f"{'Delta-stepping Time:':<25} {delta_pl_duration:.20f}")

    assert(delta_stepping_pl_result == bellman_ford_result)
    print("Correctness checks passed!")
