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
#include <map>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

void radius_stepping_serial(int source, vector<vector<edge>> &graph, vector<int> &distances, int radius) { //todo change radius to a function
    int n = graph.size();
    distances.assign(n, INF);

    distances[source] = 0;

    vector<vector<int>> S;
    S.resize(n);

    //#pragma omp for
    for (auto &e : graph[source]) {
        distances[e.to] = e.cost;
        //std::cout << "to,cost " << e.to << " " << e.cost<< "\n";
    }

    S[0].push_back(source);
    int i = 1;

    while (S[i-1].size() < n) {

        int di = INF;
        //#pragma omp for
        for (int v = 0; v < n; v++) {
            if (std::find(S[i-1].begin(), S[i-1].end(), v) == S[i-1].end()) {
                if (distances[v] != INF) {
                    //std::cout << distances[v] << "\n";
                    // #pragma omp critical
                    // {
                    di = min(di, distances[v] + radius);
                    //}
                }
                //di = min(di, distances[v] + radius);
                //std::cout << di << "\n";
            }
        }
        //std::cout << "di,INF" << di << " " << INF << "\n";

        bool repeat = true;

        while (repeat) {
            repeat = false;

            //#pragma omp for
            for (int u = 0; u < n; u++) {
                if (std::find(S[i-1].begin(), S[i-1].end(), u) == S[i-1].end() && distances[u] < di) {

                    for (auto &e : graph[u]) {
                        int v = e.to;
                        if (std::find(S[i-1].begin(), S[i-1].end(), v) == S[i-1].end()) {
                            int w = e.cost;

                            if (distances[v] > distances[u] + w) { //&& distances[u] != INF) {
                                // #pragma omp critical
                                // {
                                distances[v] = distances[u] + w;
                                //}
                                if (distances[v] <= di) {
                                    //#pragma omp critical
                                    repeat = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //#pragma omp for
        for (int v = 0; v < n; v++) {
            if (distances[v] < di) {
                // #pragma omp critical
                // {
                S[i].push_back(v);
                //}
            }
        }
        i += 1;
    }
    //return distances;
}

void radius_stepping_serial_v2_wip(int source, vector<vector<edge>> &graph, vector<int> &distances, int radius) {
    int n = graph.size();
    distances.assign(n, INF);
    distances[source] = 0;

    int i = 1;

    map<int, int> Q;
    map<int, int> R;
    for (auto &e : graph[source]) {
        Q.insert({e.cost, e.to});
        R.insert({e.cost + radius, e.to});
    } 
    Q.insert({0, 0});
    R.insert({radius, 0});

    map<int, int>::iterator itr;
    // for (itr = Q.begin(); itr != Q.end(); itr++) {
    //     cout << (*itr).first << " " << (*itr).second << "\n";
    // }

    while (Q.size() > 0) {
        auto di_pair = *(R.begin());
        int di = di_pair.first;
        //R.erase(R.begin());
        R.erase(di);

        //set<int> Ai = set<int>(Q.begin(), Q.end());
        map<int, int> Ai = map<int, int>(Q.begin(), Q.end());

        Ai.erase(Ai.lower_bound(di), Ai.end());
        Q.erase(Q.begin(), Q.lower_bound(di));

        // Ai.erase(Ai.lower_bound(di), Ai.end());
        // Q.erase(Q.begin(), Q.lower_bound(make_pair(di, 0)));

        for (auto u : Ai) {
            //R.erase(u.first);
            for (auto pair : R) {
                if (pair.second == u.second) {
                    R.erase(pair.first);
                    break;
                }
            }

            bool repeat = true;
            while (repeat) {
                repeat = false;
                for (auto pair : Ai) {
                    int u = pair.second;
                    for (auto &e : graph[u]) {
                        int v = e.to;
                        // if (u == 0) {
                        //     cout << "V: " << v << " " << u << "\n";
                        // }
                        int w = e.cost;
                        // if (distances[u] != INF) {
                        //     std::cout << distances[u] << "\n";
                        // }

                        if (distances[u] != INF && distances[v] > distances[u] + w) {
                            //std::cout << "TEST" << "\n";
                            if (distances[v] > di && distances[u] + w <= di) {
                                //R.erase(distances[v] + radius);
                                for (auto pair : R) {
                                    if (pair.second == v) {
                                        R.erase(pair.first);
                                        break;
                                    }
                                }
                                //Q.erase(distances[v]);
                                for (auto pair : Q) {
                                    if (pair.second == v) {
                                        R.erase(pair.first);
                                        break;
                                    }
                                }
                                Ai.insert({distances[v], v});
                            }
                            distances[v] = distances[u] + w;
                            // if (std::find(Ai.begin(), Ai.end(), v) != Ai.end()) {
                            //     repeat = true;
                            // }
                            for (auto pair : Ai) {
                                if (pair.second == v) {
                                    repeat = true;
                                    break;
                                }
                            }
                            // if (Ai.count(distances[v])) {
                            //     repeat = true;
                            // }


                            if (distances[v] > di) {
                                //Q.erase(distances[v]);
                                for (auto pair : Q) {
                                    if (pair.second == v) {
                                        R.erase(pair.first);
                                        break;
                                    }
                                }
                                Q.insert({distances[v], v});
                                //R.erase(distances[v] + radius);
                                for (auto pair : R) {
                                    if (pair.second == v) {
                                        R.erase(pair.first);
                                        break;
                                    }
                                }
                                R.insert({distances[v] + radius, v});
                            }
                        }
                    }  
                }
            }
        } 
        i += 1;
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

    // omp_set_num_threads(4);
    // #pragma omp parallel
    // {
    //     #pragma omp single
    //     {
    //         int num_threads = omp_get_num_threads();
    //         std::cout << "Number of threads in parallel region: " << num_threads << std::endl;
    //     }
    // }
    // printf("%d\n", correctness_check());

    if (input == "NY") {
        file_path = "data/USA-road-d.NY.gr";
    } else if (input == "CAL") {
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

    vector<int> distances_delta_p;
    int radius = 6000;

    auto start_time = std::chrono::high_resolution_clock::now();
    radius_stepping_serial(0, graph, distances_delta_p, radius);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Print the time taken
    printf("Time taken for radius_stepping_serial: %f\n", duration);

    return 0;
}
