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
#include <mutex>
#include <unordered_set>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

#define N 300000

int max_bucket;
unordered_set<int> B[N];
vector<int> S;
struct req{int w,d;};
vector<req> REQ;
std::vector<std::mutex> B_mutex(N);
std::vector<std::mutex>* distances_mutex;

vector<vector<req>> REQ_local(omp_get_max_threads());
vector<vector<int>> S_local(omp_get_max_threads());


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
    (*distances_mutex)[w].lock();
    if (d < distances[w]) {
        distances[w] = d;
        int bucket_index = distances[w] / delta;
        (*distances_mutex)[w].unlock();
        
        B_mutex[bucket_index].lock();
        B[bucket_index].erase(w);
        B_mutex[bucket_index].unlock();

        B_mutex[d / delta].lock();
        B[d / delta].insert(w);
        B_mutex[d / delta].unlock();

        if (d / delta > max_bucket) max_bucket = d / delta;
    } else {
        (*distances_mutex)[w].unlock();
    }
}


void delta_stepping_parallel(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta) {
    max_bucket = 0;
    int n = graph.size();
    distances.assign(n, INF);
    relax(source, 0, distances, delta);

    int j = 0;
    while (!bempty(j)) {
        S.clear();

        while (!B[j].empty()) {
            REQ.clear();

            #pragma omp parallel
            {
                int tid = omp_get_thread_num();

                #pragma omp single
                for (std::unordered_set<int>::iterator it = B[j].begin(); it != B[j].end(); ++it) {
                    int v = *it;
                    for (int k = 0; k < graph[v].size(); k++) {
                        if (graph[v][k].cost <= delta) {
                            req r;
                            r.w = graph[v][k].to;
                            r.d = distances[v] + graph[v][k].cost;
                            REQ_local[tid].push_back(r);
                        }
                    }
                    S_local[tid].push_back(v);
                }

                // Merge REQ_local and S_local into REQ and S using a single thread
                #pragma omp single
                {
                    for (const auto &local_reqs : REQ_local) {
                        REQ.insert(REQ.end(), local_reqs.begin(), local_reqs.end());
                    }

                    for (const auto &local_s : S_local) {
                        S.insert(S.end(), local_s.begin(), local_s.end());
                    }

                    for (auto &local_reqs : REQ_local) {
                        local_reqs.clear();
                    }

                    for (auto &local_s : S_local) {
                        local_s.clear();
                    }

                    B[j].clear();
                }

                // Relax light edges
                #pragma omp for
                for (int i = 0; i < REQ.size(); i++) {
                    relax(REQ[i].w, REQ[i].d, distances, delta);
                }
            }

            // Clear REQ
            REQ.clear();
        }

        // Heavy edges
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();

            #pragma omp for
            for (int i = 0; i < S.size(); i++) {
                int v = S[i];
                for (int k = 0; k < graph[v].size(); k++) {
                    if (graph[v][k].cost > delta) {
                        relax(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
                    }
                }
            }
        }
        j++;
    }
}

void delta_stepping_parallel_local(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
{

    max_bucket = 0;
    int n = graph.size();
    distances.assign(n, INF);
    relax(source, 0, distances, delta);
    
    int j = 0;
    while (!bempty(j)) {
        S.clear();

        while (!B[j].empty()) {
            
            #pragma omp parallel
            {
                std::vector<int> vector_j(B[j].begin(), B[j].end());
                int tid = omp_get_thread_num();
                
                #pragma omp for
                for (int i = 0; i < vector_j.size(); i++) {
                    int v = vector_j[i];
                    for (int k = 0; k < graph[v].size(); k++) {
                        if (graph[v][k].cost <= delta) {
                            req r;
                            r.w = graph[v][k].to;
                            r.d = distances[v] + graph[v][k].cost;
                            REQ_local[tid].push_back(r);
                        }
                    }
                    S_local[tid].push_back(v);
                }


                #pragma omp single
                B[j].clear();

                for (int i = 0; i < REQ_local[tid].size(); i++) {
                    relax(REQ_local[tid][i].w, REQ_local[tid][i].d, distances, delta);
                }
                
                REQ_local[tid].clear();
            }


            for (const auto &local_s : S_local) {
                S.insert(S.end(), local_s.begin(), local_s.end());
            }

            for (auto &local_s : S_local) {
                local_s.clear();
            }
        }

        // Heavy edges
        #pragma omp parallel for
        for (int i = 0; i < S.size(); i++) {
            int v = S[i];
            for (int k = 0; k < graph[v].size(); k++) {
                if (graph[v][k].cost > delta) {
                    relax(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
                }
            }
        }
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
    delta_stepping_parallel_local(0, graph, delta_stepping_distances, delta);

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

    omp_set_num_threads(4);

    // printf("%d\n", correctness_check());

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

    distances_mutex = new std::vector<std::mutex>(num_vertices);
    vector<int> distances;
    vector<int> distances_delta_s;
    vector<int> distances_delta_p;
    int delta = 200;

    // Call the delta_stepping() function
    // double t_start_p = omp_get_wtime();
    // delta_stepping_serial(0, graph, distances_delta, delta);
    // double delta_time_serial = omp_get_wtime() - t_start_p;

    // Call the dijkstra() function
    double t_dijk_start = omp_get_wtime();
    dijkstra(0, graph, distances);
    double dijks_time = omp_get_wtime() - t_dijk_start;

    double t_start_p = omp_get_wtime();
    delta_stepping_parallel_local(0, graph, distances_delta_p, delta);
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
