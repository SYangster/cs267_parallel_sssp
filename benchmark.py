from src.serial.bellman_ford import bellman_ford
from src.serial.dijkstra import dijkstra
from serial.delta_stepping import delta_stepping

import networkx as nx
import time
import random


NUM_VERTICES = 300


g = nx.complete_graph(NUM_VERTICES)
G = nx.Graph()

for u, v in g.edges:
    G.add_edge(u, v, weight=random.uniform(1, 10))


if __name__ == '__main__':
    source = 0
    delta = 5

    bf_start = time.time()
    bellman_ford_result = bellman_ford(G, source)
    bf_duration = time.time() - bf_start

    dij_start = time.time()
    dijkstra_result = dijkstra(G, source)
    dij_duration = time.time() - dij_start

    delta_step_start = time.time()
    delta_stepping_result = delta_stepping(G, source, delta)
    delta_duration = time.time() - delta_step_start

    print("Bellman-Ford Result:", bf_duration)
    print("Dijkstra Result:", dij_duration)
    print("Delta-stepping Result:", delta_duration)
