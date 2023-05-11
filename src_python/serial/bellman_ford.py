def bellman_ford(G, source):
    distance = {node: float('inf') for node in G.nodes}
    distance[source] = 0
    
    for _ in range(G.number_of_nodes() - 1):
        for u, v, data in G.edges(data=True):
            w = data['weight']
            if distance[u] + w < distance[v]:
                distance[v] = distance[u] + w
            if distance[v] + w < distance[u]:
                distance[u] = distance[v] + w
    
    for u, v, data in G.edges(data=True):
        w = data['weight']
        if distance[u] + w < distance[v] or distance[v] + w < distance[u]:
            return None # Negative cycle detected
    
    return distance