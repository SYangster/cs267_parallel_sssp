import heapq

def dijkstra(G, source):
    distance = {node: float('inf') for node in G.nodes}
    distance[source] = 0
    priority_queue = [(0, source)]
    
    while priority_queue:
        dist, u = heapq.heappop(priority_queue)
        for v, data in G[u].items():
            w = data['weight']
            if distance[u] + w < distance[v]:
                distance[v] = distance[u] + w
                heapq.heappush(priority_queue, (distance[v], v))
                
    return distance