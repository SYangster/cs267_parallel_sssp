#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include <memory>
#include <cmath> // for sqrt function
#include <chrono>
#include <omp.h>
#include <atomic>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

#define N 10000

int max_bucket;
vector<vector<int>> B;
vector<int> S;
vector<int>* vertex_to_bucket;
vector<int> bucket_sizes[N];
struct req{int w,d;};
vector<req> REQ;

//vector<vector<req>> REQ_local(omp_get_max_threads());
vector<vector<int>> S_local(omp_get_max_threads());
vector<vector<vector<int>>> B_local(omp_get_max_threads(), std::vector<std::vector<int>>(N));
vector<int> B_index;


bool bempty(int j){
	for(int i=j;i<=max_bucket;i++)
		//if(!B[i].empty())
        if (B_index[i] > 0)
			return false;
	return true;
}

// w is the vertice to go to; d is the distance
void relax(int w, int d, vector<atomic<int>> &distances, int delta, vector<vector<int>> &B_local) {
    #pragma omp critical
    {
    int d_w = distances[w].load(memory_order_relaxed);
    if (d < d_w) {
        int bucket_idx = d/delta;

        if (d_w != INF) { //&& vertex_to_bucket[w] != bucket_idx){
            
            vector<int>::iterator res = find(B_local[d_w / delta].begin(), B_local[d_w / delta].end(), w);
            if (res != B_local[d_w].end()) {
                B_local[d_w / delta].erase(res);
            }
        }
        
        //B[bucket_idx].push_back(w);
        B_local[bucket_idx].push_back(w);
        printf("bucket-id %d w: %d d: %d\n", bucket_idx, w, d);
        if (d/delta > max_bucket) max_bucket = d/delta;
        distances[w].store(d, memory_order_relaxed);
    }
    }
}

void merge_thread_buckets(int bucket_idx, bool is_light_edge, vector<vector<int>> &B_local) {
    int num_buckets = max_bucket+1;
    if (is_light_edge) {
        num_buckets = 1;
    }
    for (int i = bucket_idx; i < bucket_idx + num_buckets; i++) {
        if (!B_local[i].size()) continue;

        #pragma omp atomic
        B_index[i] += B_local[i].size();

        printf("Merging: index-%d total-%d %d\n", i, B_index[i], B_local[i][0]);
        std::copy(B_local[i].begin(), B_local[i].end(), B[i].begin() + B_index[i] - B_local[i].size());
    }
}

void delta_stepping_parallel(int source, vector<vector<edge>> &graph, vector<atomic<int>> &distances, int delta)
{   

    max_bucket = 0;
    distances[source].store(0, memory_order_relaxed);
    B[0].push_back(source);
    B_index[source] = 1;
    

    int j = 0;
    while(!bempty(j)){
        //S.clear();
        // printf("%d \n", j);
        //while (!B[j].empty())
        while (B_index[j] != 0) {

            printf("running index: %d \n", j);

            #pragma omp parallel 
            {
                
                int tid = omp_get_thread_num();
                // look for vertices within delta
                //for (int i = 0; i < B[j].size(); i++){
                #pragma omp for
                for (int i = 0; i < B_index[j]; i++){
                    int v = B[j][i];
                    for (int k=0; k<graph[v].size(); k++){
                        // add requests: from v to light neighbors
                        if (graph[v][k].cost <= delta){
                            req r;
                            r.w = graph[v][k].to;
                            r.d = distances[v] + graph[v][k].cost;
                            printf("light %d %d \n", r.w, r.d);
                            printf("%d %d \n", distances[r.w].load());
                            relax(r.w, r.d, distances, delta, B_local[tid]);
                        }
                    }
                    S_local[tid].push_back(v);
                }

                #pragma omp single
                B_index[j] = 0;
                //B[j].clear();

                merge_thread_buckets(j, true, B_local[tid]);

                if (j == 1) {
                    printf("B-indx %d \n", B_index[1]);
                    for (int i = 0; i < B_index[1]; i++){
                        printf("B: %d ", B[j][i]);
                    }
                }


                #pragma omp single
                B_local[tid][j].clear();
                
            }
            
        }
        

        // if (j == 1) {
        //     for (int i = 0; i < B_index[2]; i++){
        //         printf("here %d %d \n", B[2][i], j);
        //     }
        //     }

        #pragma omp parallel 
        {
            int tid = omp_get_thread_num();

            // printf("%d \n", S_local[tid].size());
            // heavy edge
            for (int i = 0; i < S_local[tid].size(); i++) {
                int v = S_local[tid][i];
                for (int k=0; k<graph[v].size(); k++){
                    // heavy edge requests
                    if (graph[v][k].cost > delta) {
                        req r;
                        r.w = graph[v][k].to;
                        r.d = distances[v] + graph[v][k].cost;
                        printf("S_heavy %d %d \n", r.w, r.d);
                        relax(r.w, r.d, distances, delta, B_local[tid]);
                    }
                }
            }
            // if (j == 1) {printf("2222 %d \n", B_local[0][2][0]);}
            merge_thread_buckets(j, false, B_local[tid]);
            // if (j == 1) {printf("3333 %d \n", B[2][0]);}
            S_local[tid].clear();
        }

        // for (int i = 0; i < 10; i++)
        // {
        //     printf("%d ", B_index[i]);
        // }
        // if (true) {
        //     for (int i = 0; i < B_index[j+1]; i++){
        //         printf("Final Bucket at index: %d %d \n", j, B[j+1][i]);
        //     }
        //     }
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

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset>" << endl;
        return 1;
    }

    string input = argv[1];
    string file_path;

    // printf("%d\n", correctness_check());
    omp_set_num_threads(1);

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
    
    B_index.resize(N, 0);
    B.resize(N);
    for(int i = 0; i < N; ++i) {
        B[i].resize(num_vertices, 0);
    }

    vertex_to_bucket = new vector<int>(num_vertices);


    vector<atomic<int>> distances_delta_p(num_vertices);
    for(auto& v : distances_delta_p) {
        v.store(INF, std::memory_order_relaxed);
    }
    int delta = 2000;

    vector<int> distances;
    double t_dijk_start = omp_get_wtime();
    dijkstra(0, graph, distances);
    double dijks_time = omp_get_wtime() - t_dijk_start;

    double t_start_p = omp_get_wtime();
    delta_stepping_parallel(0, graph, distances_delta_p, delta);
    double delta_time_parallel = omp_get_wtime() - t_start_p;


    for (size_t i = 0; i < distances.size(); ++i)
    {
        if (distances[i] != distances_delta_p[i].load())
        {   
            cout << i  << "\t" << distances[i] << "\t" << distances_delta_p[i] << "\n";
            return 0;
        }
    }

    // auto start_time = std::chrono::high_resolution_clock::now();
    // delta_stepping_serial(0, graph, distances_delta_p, delta);
    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Print the time taken
    printf("Time taken for delta_stepping serial: %f\n", delta_time_parallel);

    return 0;
}
