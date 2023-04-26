#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include <memory>
#include <cmath> // for sqrt function
#include <tuple>
#include <chrono>
#include "load_graph.h"


const int INF = numeric_limits<int>::max();

void dijkstra(int source, vector<vector<edge>> &graph, vector<int> &distances)
{
    int n = graph.size();
    distances.assign(n, INF);
    distances[source] = 0;

    // Use a set to keep track of unprocessed vertices
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> unprocessed_vertices;
    unprocessed_vertices.push({0, source});

    while (!unprocessed_vertices.empty())
    {
        // Get the vertex with the smallest distance from the set
        auto u = unprocessed_vertices.top();
        unprocessed_vertices.pop();

        int vertex_u = u.second;

        if (distances[vertex_u] < u.first) continue;
        // Relax all the outgoing edges from the vertex
        for (auto &e : graph[vertex_u])
        {
            int v = e.to;
            int weight = e.cost;
            int new_distance = distances[vertex_u] + weight;
            if (new_distance < distances[v])
            {
                distances[v] = new_distance;
                unprocessed_vertices.push({distances[v], v});
            }
        }
    }
}


int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset>" << endl;
        return 1;
    }

    string input = argv[1];
    string file_path;

    // printf("%d\n", correctness_check());

    if (input == "NY") {
        file_path = "data/USA-road-d.NY.gr";
    } else if (input == "CAL")
    {
        file_path = "data/USA-road-d.CAL.gr";
    } else if (input == "W") {
        file_path = "data/USA-road-d.W.gr";
    } else if (input == "USA") {
        file_path = "data/USA-road-d.USA.gr";
    }
    
    string dataset_name;
    int num_vertices, num_edges;
    vector<vector<edge>> graph;

    tie(dataset_name, num_vertices, num_edges, graph) = load_gr_file(file_path);

    cout << "Dataset: " << file_path << endl;
    cout << "Number of vertices: " << num_vertices << endl;
    cout << "Number of edges: " << num_edges << endl;

    vector<int> distances;

    // Call the dijkstra() function
    auto start_time = std::chrono::high_resolution_clock::now();
    dijkstra(0, graph, distances);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Print the time taken
    printf("Time taken for dijkstra: %f\n", duration);

    return 0;
}
