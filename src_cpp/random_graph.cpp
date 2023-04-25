#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>

using namespace std;
struct edge {
    int to, cost;
};

vector<vector<edge>> generate_random_graph(int num_vertices, double edge_density, int min_weight, int max_weight) {
    srand(time(0));
    vector<vector<edge>> graph(num_vertices);

    int max_edges = num_vertices * (num_vertices - 1) / 2;
    int num_edges = static_cast<int>(edge_density * max_edges);

    for (int i = 0; i < num_edges; ++i) {
        int u = rand() % num_vertices;
        int v = rand() % num_vertices;
        int weight = rand() % (max_weight - min_weight + 1) + min_weight;

        // Ensure the edge doesn't already exist and u != v
        auto it = find_if(graph[u].begin(), graph[u].end(), [v](const edge& e) { return e.to == v; });
        if (u != v && it == graph[u].end()) {
            graph[u].push_back({v, weight});
            graph[v].push_back({u, weight});
        } else {
            --i; // Decrement the counter to try again
        }
    }

    return graph;
}

