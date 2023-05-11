
import random
import math
import networkx as nx
from collections import defaultdict, deque
from multiprocessing import Pool, Manager, Lock

def parallel_relax(args):
    u, light, heavy, delta, distance, buckets = args
    local_updates = []
    for v, w in light[u]:
        relaxed = relax(u, v, w, delta, distance, buckets)
        if relaxed:
            local_updates.append((v, distance[v]))
    for v, w in heavy[u]:
        relaxed = relax(u, v, w, delta, distance, buckets)
        if relaxed:
            local_updates.append((v, distance[v]))
    return local_updates

def relax(u, v, w, delta, distance, buckets):
    if distance[u] + w < distance[v]:
        distance[v] = distance[u] + w
        return True
    return False

def delta_stepping_multiprocessing(G, source, delta):
    n = G.number_of_nodes()
    max_distance = n * max(data['weight'] for u, v, data in G.edges(data=True))
    num_buckets = math.ceil(max_distance / delta)
    heavy = defaultdict(list)
    light = defaultdict(list)
    
    for u, v, data in G.edges(data=True):
        w = data['weight']
        if w > delta:
            heavy[u].append((v, w))
            heavy[v].append((u, w))
        else:
            light[u].append((v, w))
            light[v].append((u, w))

    distance = [float('inf')] * n
    distance[source] = 0
    buckets = [deque() for _ in range(num_buckets)]
    buckets[0].append(source)
    
    lock = Lock()
    
    for r in range(num_buckets):
        while buckets[r]:
            with Pool() as pool:
                results = pool.map(parallel_relax, [(u, light, heavy, delta, distance.copy(), buckets) for u in buckets[r]])
                buckets[r].clear()
                with lock:
                    for local_updates in results:
                        for v, dist in local_updates:
                            if dist < distance[v]:
                                distance[v] = dist
                                bucket_idx = int(distance[v] // delta)
                                buckets[bucket_idx].append(v)
        
    return {node: distance[node] for node in range(n)}