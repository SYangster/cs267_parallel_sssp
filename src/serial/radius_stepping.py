from collections import defaultdict, deque
import networkx as nx

# def radius_stepping_precompute(G):
#     r = {node: float('inf') for node in G.nodes}
#     p = 3
#     for u in G.nodes():
#         distance_to_pth_closest = float('inf')
#         for _ in range(p):
#     return r

def radius_stepping(G, source, r):
    n = G.number_of_nodes()
    distance = {node: float('inf') for node in G.nodes}
    distance[source] = 0

    S = {node: [] for node in G.nodes}

    for v in nx.neighbors(G, source):
        distance[v] = G.get_edge_data(source, v)["weight"]

    S[0] = [source]
    i = 1

    while len(S[i-1]) < n:
        di = min (map(lambda v : distance[v] + r, [v for v in G.nodes() if v not in S[i-1]]))

        repeat = True
        while repeat:
            repeat = False
            for u in [u for u in G.nodes() if u not in S[i-1] and distance[u] < di]:
                for v in [v for v in nx.neighbors(G, u) if v not in S[i-1]]:
                    w = G.get_edge_data(u, v)["weight"]

                    if distance[v] > distance[u] + w:
                        distance[v] = distance[u] + w
                        if distance[v] <= di:
                            repeat = True

        S[i] = [v for v in G.nodes if distance[v] < di]
        i += 1
    return distance
