// random_graph.h
#pragma once

#include <vector>
#include <utility>

using namespace std;
struct edge {
    int to, cost;
};

vector<vector<edge>> generate_random_graph(int num_vertices, double edge_density, int min_weight, int max_weight);
