#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <tuple>

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

tuple<string, int, int, vector<vector<edge>>> load_gr_file(const string& file_path) {
    ifstream file(file_path);
    if (!file.is_open()) {
        cerr << "Error opening file: " << file_path << endl;
        exit(1);
    }

    string line;
    string dataset_name;
    int num_vertices, num_edges;

    // Read header lines and find the line with problem descriptor
    while (getline(file, line)) {
        if (line[0] == 'c') {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> dataset_name;
        }
        if (line[0] == 'p') {
            istringstream iss(line);
            string tmp;
            iss >> tmp >> tmp >> num_vertices >> num_edges;
            break;
        }
    }

    vector<vector<edge>> graph(num_vertices);

    // Read edges
    while (getline(file, line)) {
        if (line[0] == 'a') {
            istringstream iss(line);
            char tmp;
            int from, to, cost;
            iss >> tmp >> from >> to >> cost;

            // Indices in .gr files are 1-based, we need 0-based indices
            from--;
            to--;

            graph[from].push_back({to, cost});
            // Assuming the graph is undirected
            graph[to].push_back({from, cost});
        }
    }

    return make_tuple(dataset_name, num_vertices, num_edges, graph);
}

