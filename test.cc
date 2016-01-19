#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <random>
#include "adaptive-betweenness-centrality.h"
#include <gflags/gflags.h>
using namespace std;
using namespace adaptive_betweenness_centrality;

DEFINE_string(G, "", "input file");
DEFINE_string(method, "coverage", "coverage / betweenness");
DEFINE_int32(M, 1024, "M");
DEFINE_int32(k, 2, "k");


int main(int argc, char *argv[])
{
  srandom(time(NULL));
  google::ParseCommandLineFlags(&argc, &argv, true);

  cout << "# G = " << FLAGS_G << endl;
  cout << "# method = " << FLAGS_method << endl;
  cout << "# M = " << FLAGS_M << endl;
  cout << "# k = " << FLAGS_k << endl;

  g_t G;
  load_graph(FLAGS_G, G);
  int V = G.size();

  cout << "# V = " << V << endl;

  vector<int> seeds;
  if (FLAGS_method == "exact-coverage") {
    vector<int> btwss;
    exact_coverage(G, seeds, btwss);
    for (int i = 0; i < V; i++) {
      cout << i << " " << btwss[i] << endl;
    }
  } else if (FLAGS_method == "approximate-coverage") {
    vector<int> btwss;
    approximate_coverage(G, FLAGS_M, seeds, btwss);
    for (int i = 0; i < V; i++) {
      cout << i << " " << btwss[i] << endl;
    }
  } else if (FLAGS_method == "topk-coverage") {
    vector<int> btwss;
    adaptive_approximate_coverage(G, FLAGS_M, FLAGS_k, seeds, btwss);
    for (int i = 0; i < (int)seeds.size(); i++) {
      cout << i << " " << btwss[i] << endl;
    }
  } else if (FLAGS_method == "exact-betweenness") {
    vector<double> btwss;
    exact_betweenness(G, seeds, btwss);
    for (int i = 0; i < V; i++) {
      cout << i << " " << btwss[i] << endl;
    }
  } else if (FLAGS_method == "approximate-betweenness") {
    vector<double> btwss;
    approximate_betweenness(G, FLAGS_M, seeds, btwss);
    for (int i = 0; i < V; i++) {
      cout << i << " " << btwss[i] << endl;
    }
  } else if (FLAGS_method == "topk-betweenness") {
    vector<double> btwss;
    adaptive_approximate_betweenness(G, FLAGS_M, FLAGS_k, seeds, btwss);
    for (int i = 0; i < (int)seeds.size(); i++) {
      cout << i << " " << btwss[i] << endl;
    }
  }
}