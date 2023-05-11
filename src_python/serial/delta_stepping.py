from collections import defaultdict, deque

def delta_stepping(G, source, delta):
    n = G.number_of_nodes()
    distance = {node: float('inf') for node in G.nodes}
    distance[source] = 0
    heavy = defaultdict(list)
    light = defaultdict(list)
    buckets = [deque() for _ in range(n)]
    
    for u, v, data in G.edges(data=True):
        w = data['weight']
        if w > delta:
            heavy[u].append((v, w))
            heavy[v].append((u, w))
        else:
            light[u].append((v, w))
            light[v].append((u, w))
    
    def relax(u, v, w, bucket_idx):
        if distance[u] + w < distance[v]:
            distance[v] = distance[u] + w
            buckets[bucket_idx].append(v)
    
    buckets[0].append(source)
    for r in range(n):
        while buckets[r]:
            u = buckets[r].popleft()
            for v, w in light[u]:
                relax(u, v, w, r)
            for v, w in heavy[u]:
                relax(u, v, w, int((distance[u] + w) // delta))
    
    return distance