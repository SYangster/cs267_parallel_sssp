#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <set>
#include <memory>
#include <cmath> // for sqrt function
#include <chrono>
#include "load_graph.h"

using namespace std;
const int INF = numeric_limits<int>::max();

#define N 1408199

int max_bucket;
vector<int> B[N];
vector<int> S;
struct req{int w,d;};
vector<req> REQ;


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

void delta_stepping_serial(int source, vector<vector<edge>> &graph, vector<int> &distances, int delta)
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
            for (int i=0; i<REQ.size(); i++){
                relax(REQ[i].w, REQ[i].d, distances, delta);
            }

        }
        
        REQ.clear();

        // heavy edge
        for (int i=0; i<S.size(); i++) {
            int v = S[i];
            for (int k=0; k<graph[v].size(); k++){
                // heavy edge requests
                if (graph[v][k].cost > delta) {
                    req r;
                    r.w = graph[v][k].to;
                    r.d = distances[v] + graph[v][k].cost;
                    if (j == 1) printf("S_heavy %d %d \n", r.w, r.d);
                    REQ.push_back(r);
                }
            }
        }
     
        for (int i=0; i<REQ.size(); i++){
            relax(REQ[i].w, REQ[i].d, distances, delta);
        }
        if (j == 1) {
                for (int i: B[2]){
                    printf("here %d %d \n", i, j);
                }
            }
        j++;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <dataset>" << endl;
        return 1;
    }

    string input = argv[1];
    string input2 = argv[2];
    int delta = stoi(input2);
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

    vector<int> distances_delta_p;
    // int delta = 2000;

    auto start_time = std::chrono::high_resolution_clock::now();
    delta_stepping_serial(0, graph, distances_delta_p, delta);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Print the time taken
    printf("Time taken for delta_stepping serial: %f\n", duration);

    return 0;
}
