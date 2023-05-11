#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include <memory>
#include <cmath> // for sqrt function
#include <omp.h>
#include <chrono>
#include <mutex>
#include <unordered_set>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

#define N 1408199

int max_bucket;
deque<int> B[N];
vector<int> S;
struct req{int w,d;};
vector<req> REQ;
std::vector<std::mutex> B_mutex(N);
std::vector<std::mutex>* distances_mutex;

vector<vector<req>> REQ_local(omp_get_max_threads());
vector<vector<int>> S_local(omp_get_max_threads());


bool bempty(int j){
	for(int i=j;i<=max_bucket;i++)
		if(!B[i].empty())
			return false;
	return true;
}

void relax(int w, int d, vector<int> &distances, int delta) {
    #pragma omp critical
    {
        if (d < distances[w]) {
            if (distances[w] != INF) {
                B[distances[w] / delta].erase(remove(B[distances[w] / delta].begin(), B[distances[w] / delta].end(), w), B[distances[w] / delta].end());
            }
            B[d / delta].push_back(w);
            if (d / delta > max_bucket) max_bucket = d / delta;
            distances[w] = d;
        }
    }
}

// void relax_lock(int w, int d, vector<int> &distances, int delta) {
//     (*distances_mutex)[w].lock();
//     if (d < distances[w]) {
//         distances[w] = d;
//         int bucket_index = distances[w] / delta;
//         (*distances_mutex)[w].unlock();
        
//         B_mutex[bucket_index].lock();
//         B[bucket_index].erase(w);
//         B_mutex[bucket_index].unlock();

//         B_mutex[d / delta].lock();
//         B[d / delta].insert(w);
//         B_mutex[d / delta].unlock();

//         if (d / delta > max_bucket) max_bucket = d / delta;
//     } else {
//         (*distances_mutex)[w].unlock();
//     }
// }


// void delta_stepping_parallel(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
// {
//     max_bucket = 0;
//     int n = graph.size();
//     distances.assign(n, INF);
//     relax(source, 0, distances, delta);
    
//     int j = 0;
//     while(!bempty(j)){
//         S.clear();

//         while (!B[j].empty())
//         {
//             REQ.clear();

//             // look for vertices within delta
//             #pragma omp parallel for
//             for (int i = 0; i < B[j].size(); i++){
//                 int v = B[j][i];
//                 for (int k=0; k<graph[v].size(); k++){
//                     // add requests: from v to light neighbors
//                     if (graph[v][k].cost <= delta){
//                         req r;
//                         r.w = graph[v][k].to;
//                         r.d = distances[v] + graph[v][k].cost;
//                         #pragma omp critical
//                         {
//                             REQ.push_back(r);
//                         }
                        
//                     }
//                 }
//                 #pragma omp critical
//                 {
//                     S.push_back(v);
//                 }
//             }

//             B[j].clear();
            
//             // relax light edges
//             #pragma omp parallel for
//             for (int i=0; i<REQ.size(); i++){
//                 relax(REQ[i].w, REQ[i].d, distances, delta);
//             }

//         }
//         REQ.clear();

//         // heavy edges
//         #pragma omp parallel for
//         for (int i = 0; i < S.size(); i++) {
//             int v = S[i];
//             for (int k = 0; k < graph[v].size(); k++) {
//                 if (graph[v][k].cost > delta) {
//                     relax(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
//                 }
//             }
//         }
//         j++;
//     }
// }

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
                int tid = omp_get_thread_num();

                #pragma omp for
                for (int i = 0; i < B[j].size(); i++) {
                    int v = B[j][i];
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
            }

            B[j].clear();

            for (int t = 0; t < omp_get_max_threads(); t++) {
                for (int i = 0; i < REQ_local[t].size(); i++) {
                    relax(REQ_local[t][i].w, REQ_local[t][i].d, distances, delta);
                }
                REQ_local[t].clear();
            }
        }


        for (auto &local_s : S_local) {
            for (int i = 0; i < local_s.size(); i++){
                int v = local_s[i];
                for (int k = 0; k < graph[v].size(); k++) {
                    if (graph[v][k].cost > delta) {
                        relax(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
                    }
                }
            }
            local_s.clear();
        }
        // for (const auto &local_s : S_local) {
        //     S.insert(S.end(), local_s.begin(), local_s.end());
        // }

        // for (auto &local_s : S_local) {
        //     local_s.clear();
        // }

        // Heavy edges
        // for (int i = 0; i < S.size(); i++) {
        //     int v = S[i];
        //     for (int k = 0; k < graph[v].size(); k++) {
        //         if (graph[v][k].cost > delta) {
        //             relax(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
        //         }
        //     }
        // }
        j++;
    }
}

// void delta_stepping_parallel_locks(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
// {

//     max_bucket = 0;
//     int n = graph.size();
//     distances.assign(n, INF);
//     relax_lock(source, 0, distances, delta);
    
//     int j = 0;
//     while (!bempty(j)) {
//         S.clear();

//         while (!B[j].empty()) {
            
//             #pragma omp parallel
//             {
//                 std::vector<int> vector_j(B[j].begin(), B[j].end());
//                 int tid = omp_get_thread_num();
                
//                 #pragma omp for
//                 for (int i = 0; i < vector_j.size(); i++) {
//                     int v = vector_j[i];
//                     for (int k = 0; k < graph[v].size(); k++) {
//                         if (graph[v][k].cost <= delta) {
//                             req r;
//                             r.w = graph[v][k].to;
//                             r.d = distances[v] + graph[v][k].cost;
//                             REQ_local[tid].push_back(r);
//                         }
//                     }
//                     S_local[tid].push_back(v);
//                 }


//                 #pragma omp single
//                 B[j].clear();

//                 for (int i = 0; i < REQ_local[tid].size(); i++) {
//                     relax_lock(REQ_local[tid][i].w, REQ_local[tid][i].d, distances, delta);
//                 }
                
//                 REQ_local[tid].clear();
//             }


//             for (const auto &local_s : S_local) {
//                 S.insert(S.end(), local_s.begin(), local_s.end());
//             }

//             for (auto &local_s : S_local) {
//                 local_s.clear();
//             }
//         }

//         // Heavy edges
//         #pragma omp parallel for
//         for (int i = 0; i < S.size(); i++) {
//             int v = S[i];
//             for (int k = 0; k < graph[v].size(); k++) {
//                 if (graph[v][k].cost > delta) {
//                     relax_lock(graph[v][k].to, distances[v] + graph[v][k].cost, distances, delta);
//                 }
//             }
//         }
//         j++;
//     }
// }





// void delta_stepping_serial(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
// {
//     max_bucket = 0;
//     int n = graph.size();
//     distances.assign(n, INF);
    
//     relax(source, 0, distances, delta);
    
//     int j = 0;
//     while(!bempty(j)){
//         S.clear();

//         while (!B[j].empty())
//         {
//             REQ.clear();

//             // look for vertices within delta
//             for (int i = 0; i < B[j].size(); i++){
//                 int v = B[j][i];
//                 for (int k=0; k<graph[v].size(); k++){
//                     // add requests: from v to light neighbors
//                     if (graph[v][k].cost <= delta){
//                         req r;
//                         r.w = graph[v][k].to;
//                         r.d = distances[v] + graph[v][k].cost;
//                         REQ.push_back(r);
//                     }
//                 }
//                 S.push_back(v);
//             }
//             B[j].clear();
            
//             // relax light edges
//             for (int i=0; i<REQ.size(); i++){
//                 relax(REQ[i].w, REQ[i].d, distances, delta);
//             }
//         }
//         REQ.clear();

//         // heavy edge
//         for (int i=0; i<S.size(); i++) {
//             int v = S[i];
//             for (int k=0; k<graph[v].size(); k++){
//                 // heavy edge requests
//                 if (graph[v][k].cost > delta) {
//                     req r;
//                     r.w = graph[v][k].to;
//                     r.d = distances[v] + graph[v][k].cost;
//                     REQ.push_back(r);
//                 }
//             }
//         }
     
//         for (int i=0; i<REQ.size(); i++)
//             relax(REQ[i].w, REQ[i].d, distances, delta);
//         j++;
//     }
// }

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset>" << endl;
        return 1;
    }

    string input = argv[1];
    string file_path;

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

    distances_mutex = new std::vector<std::mutex>(num_vertices);
    vector<int> distances_delta_p;
    int delta = 2000;

    // auto start_time1 = std::chrono::high_resolution_clock::now();
    // delta_stepping_parallel(0, graph, distances_delta_p, delta);
    // auto end_time1 = std::chrono::high_resolution_clock::now();
    // auto duration1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time1 - start_time1).count();


    // auto start_time2 = std::chrono::high_resolution_clock::now();
    // delta_stepping_serial(0, graph, distances_delta_p, delta);
    // auto end_time2 = std::chrono::high_resolution_clock::now();
    // auto duration2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time2 - start_time2).count();

    auto start_time = std::chrono::high_resolution_clock::now();
    delta_stepping_parallel_local(0, graph, distances_delta_p, delta);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Print the time taken
    // printf("Time taken for serial delta_stepping parallel: %f\n", duration2);
    // printf("Time taken for og delta_stepping parallel: %f\n", duration1);
    printf("Time taken for delta_stepping parallel: %f\n", duration);
    
    return 0;
}
