#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include <memory>
#include <cmath> // for sqrt function
#include <omp.h>
#include <tuple>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

#define N 100000

int max_bucket;
vector<int> B[N];
vector<int> S;
struct req{int w,d;};
vector<req> REQ;


#define clear_distances() { \
    for (int i = 0; i < distances.size(); i++) { \
        distances[i] = INF; \
    } \
} \

bool bempty(int j){
	for(int i=j;i<=max_bucket;i++)
		if(!B[i].empty())
			return false;
	return true;
}

// w is the vertice to go to; d is the distance
void relax(int w, int d, vector<int> &distances, int delta) {
    if (d < distances[w]) {
        if (distances[w] != INF){
            
            vector<int>::iterator res = find(B[distances[w] / delta].begin(), B[distances[w] / delta].end(), w);
            if (res != B[distances[w]/delta].end()) {
                B[distances[w] / delta].erase(res);
            }
        }

        B[d/delta].push_back(w);
        if (d/delta > max_bucket) max_bucket = d/delta;
        distances[w] = d;
    }
}

void delta_stepping_parallel(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
{   

    max_bucket = 0;
    int n = graph.size();
    distances.assign(n, INF);
    
    relax(source, 0, distances, delta);
    
    int j = 0;
    while(!bempty(j)){
        S.clear();

        while (!B[j].empty())
        {
            REQ.clear();

            // look for vertices within delta
            #pragma omp for
            for (int i = 0; i < B[j].size(); i++){
                int v = B[j][i];
                for (int k=0; k<graph[v].size(); k++){
                    // add requests: from v to light neighbors
                    if (graph[v][k].cost <= delta){
                        req r;
                        r.w = graph[v][k].to;
                        r.d = distances[v] + graph[v][k].cost;
                        REQ.push_back(r);
                    }
                }
                S.push_back(v);
            }
            B[j].clear();
            
            // relax light edges
            #pragma omp for
            for (int i=0; i<REQ.size(); i++){
                relax(REQ[i].w, REQ[i].d, distances, delta);
            }
        }
        REQ.clear();

        // heavy edges
        #pragma omp for
        for (int i=0; i<S.size(); i++) {
            int v = S[i];
            for (int k=0; k<graph[v].size(); k++){
                // heavy edge requests
                if (graph[v][k].cost > delta) {
                    req r;
                    r.w = graph[v][k].to;
                    r.d = distances[v] + graph[v][k].cost;
                    REQ.push_back(r);
                }
            }
        }
        
        #pragma omp for 
        for (int i=0; i<REQ.size(); i++)
            relax(REQ[i].w, REQ[i].d, distances, delta);
        j++;
    }
}


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


bool correctness_check()
{
    int num_vertices = 500;
    double edge_density = 0.3;
    int min_weight = 1;
    int max_weight = 100;

    // Generate a random graph
    auto graph = generate_random_graph(num_vertices, edge_density, min_weight, max_weight);

    // Run Dijkstra's algorithm
    std::vector<int> dijkstra_distances;
    dijkstra(0, graph, dijkstra_distances);

    // cout << "Vertex\tDistance from Source\n";
    // for (int i = 0; i < dijkstra_distances.size(); ++i)
    //         cout << i << "\t" << dijkstra_distances[i] << "\n";

    // Run Delta-stepping algorithm
    int delta = 50;
    std::vector<int> delta_stepping_distances;
    delta_stepping_parallel(0, graph, delta_stepping_distances, delta);

    // cout << "Vertex\tDistance from Source\n";
    // for (int i = 0; i < delta_stepping_distances.size(); ++i)
    //         cout << i << "\t" << delta_stepping_distances[i] << "\n";
    

    if (dijkstra_distances.size() != delta_stepping_distances.size())
    {
        return false;
    }

    for (size_t i = 0; i < dijkstra_distances.size(); ++i)
    {
        if (dijkstra_distances[i] != delta_stepping_distances[i])
        {
            return false;
        }
    }

    return true;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset>" << endl;
        return 1;
    }

    string input = argv[1];
    string file_path;

    omp_set_num_threads(8);

    printf("%d\n", correctness_check());

    if (input == "NY") {
        file_path = "data/USA-road-d.NY.gr";
    } else if (input == "CAL")
    {
        file_path = "data/USA-road-d.CAL.gr";
    } else if (input == "W") {
        file_path = "data/USA-road-d.W.gr";
    }
    
    string dataset_name;
    int num_vertices, num_edges;
    vector<vector<edge>> graph;

    tie(dataset_name, num_vertices, num_edges, graph) = load_gr_file(file_path);

    cout << "Dataset: " << file_path << endl;
    cout << "Number of vertices: " << num_vertices << endl;
    cout << "Number of edges: " << num_edges << endl;

    vector<int> distances;
    vector<int> distances_delta_s;
    vector<int> distances_delta_p;
    int delta = 2000;

    // Call the delta_stepping() function
    // double t_start_p = omp_get_wtime();
    // delta_stepping_serial(0, graph, distances_delta, delta);
    // double delta_time_serial = omp_get_wtime() - t_start_p;

    // Call the dijkstra() function
    double t_dijk_start = omp_get_wtime();
    dijkstra(0, graph, distances);
    double dijks_time = omp_get_wtime() - t_dijk_start;

    double t_start_p = omp_get_wtime();
    delta_stepping_parallel(0, graph, distances_delta_p, delta);
    double delta_time_parallel = omp_get_wtime() - t_start_p;


    for (size_t i = 0; i < distances.size(); ++i)
    {
        if (distances[i] != distances_delta_p[i])
        {   
            cout << i  << "\t" << distances[i] << "\t" << distances_delta_p[i] << "\n";
            return 0;
        }
    }

    // Print the time taken
    printf("Time taken for dijkstra: %f\n", dijks_time);
    printf("Time taken for delta_stepping parallel: %f\n", delta_time_parallel);

    return 0;
}
