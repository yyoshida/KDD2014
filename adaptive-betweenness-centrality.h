#pragma once

#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdarg>

#define rep(i, n) for (int i = 0; i < (int)(n); i++)

namespace adaptive_betweenness_centrality {

  const double eps = 1e-8;

  typedef std::vector<std::vector<int> > g_t;
  typedef std::vector<int> he_t;
  typedef std::vector<he_t> hg_t;
  typedef std::vector<std::pair<int, double> > whe_t;
  typedef std::vector<whe_t> whg_t;

  void load_graph(const std::string& fn, g_t& G);
  void bfs(g_t& G, int s, std::vector<int>& dists, std::vector<int>& nums);
  void bfs(g_t& G, int s, std::unordered_map<int, int>& dists, std::unordered_map<int, int>& nums, std::unordered_set<int>& alive);


  void exact_betweenness(g_t& G, std::vector<int>& seeds, std::vector<double>& btws);
  void approximate_betweenness(g_t& G, int M, std::vector<int>& seeds, std::vector<double>& btws);
  void adaptive_approximate_betweenness(g_t& G, int M, int k, std::vector<int>& seeds, std::vector<double>& btwss);

  void exact_coverage(g_t& G, std::vector<int>& seeds, std::vector<int>& btws);
  void approximate_coverage(g_t& G, int M, std::vector<int>& seeds, std::vector<int>& btws);
  void adaptive_approximate_coverage(g_t& G, int M, int k, std::vector<int>& seeds, std::vector<int>& btwss);

  void build_betweenness_hypergraph(g_t& G, whg_t& H, int M, std::vector<int>& seeds, std::vector<std::pair<int, int> >* pairs = NULL);
  void build_coverage_hypergraph(g_t& G, hg_t& H, int M, std::vector<int>& seeds);

  /*********/
  /*  BFS  */
  /*********/

  void bfs(g_t& G, int s, std::vector<int>& dists, std::vector<int>& nums) {
    int V = G.size();
    int qh = 0, qt = 0;
    std::vector<int> q(V);

    dists = std::vector<int>(V, -1);
    nums = std::vector<int>(V, 0);

    q[qt++] = s;
    dists[s] = 0;
    nums[s] = 1;
    while (qh != qt) {
      int u = q[qh++];
      rep (i, G[u].size()) {
        int v = G[u][i];

        if (dists[v] == -1 || dists[v] == dists[u] + 1) {
          nums[v] += nums[u];
        }
        if (dists[v] == -1) {
          dists[v] = dists[u] + 1;
          q[qt++] = v;
        }
      }
    }
  }

  void bfs(g_t& G, int s, std::unordered_map<int, int>& dists, std::unordered_map<int, int>& nums, std::unordered_set<int>& alive) {
    int qh = 0;
    std::vector<int> q;

    dists.clear();
    nums.clear();

    if (alive.find(s) == alive.end()) return;

    q.push_back(s);
    dists[s] = 0;
    nums[s] = 1;
    while (qh != (int)q.size()) {
      int u = q[qh++];
      rep (i, G[u].size()) {
        int v = G[u][i];
        if (alive.find(v) == alive.end()) continue;

        if (dists.find(v) == dists.end() || dists[v] == dists[u] + 1) {
          nums[v] += nums[u];
        }
        if (dists.find(v) == dists.end()) {
          dists[v] = dists[u] + 1;
          q.push_back(v);
        }
      }
    }
  }

  void restricted_bfs(g_t& G, int s, std::unordered_map<int, int>& baseline_dists, std::unordered_map<int, int>& nums, std::vector<bool>& is_seed, std::unordered_set<int>& domain) {
    int qh = 0;
    std::vector<int> q;
    std::unordered_map<int, int> dists;
    nums.clear();

    if (baseline_dists.find(s) != baseline_dists.end() && baseline_dists[s] != 0) return;
    if (domain.count(s) == 0) return;

    q.push_back(s);
    dists[s] = 0;
    nums[s] = 1;
    while (qh != (int)q.size()) {
      int u = q[qh++];
      rep (i, G[u].size()) {
        int v = G[u][i];
        if (domain.count(v) == 0 || is_seed[v]) continue;
        if (baseline_dists[v] != baseline_dists[u] + 1) continue;

        if (dists.find(v) == dists.end() || dists[v] == dists[u] + 1) {
          nums[v] += nums[u];
        }
        if (dists.find(v) == dists.end()) {
          dists[v] = dists[u] + 1;
          q.push_back(v);
        }
      }
    }
  }



  void restricted_bfs(g_t& G, int s, const std::vector<int>& baseline_dists, std::vector<int>& nums, std::vector<bool>& is_seed) {
    int V = G.size();
    int qh = 0, qt = 0;
    std::vector<int> q(V);

    std::vector<int> dists(V, -1);
    nums = std::vector<int>(V, 0);

    if (baseline_dists[s] != 0) return;

    q[qt++] = s;
    dists[s] = 0;
    nums[s] = 1;
    while (qh != qt) {
      int u = q[qh++];
      rep (i, G[u].size()) {
        int v = G[u][i];
        if (is_seed[v]) continue;
        if (baseline_dists[v] != baseline_dists[u] + 1) continue;

        if (dists[v] == -1 || dists[v] == dists[u] + 1) {
          nums[v] += nums[u];
        }
        if (dists[v] == -1) {
          dists[v] = dists[u] + 1;
          q[qt++] = v;
        }
      }
    }
  }



  /*************************/
  /*  Betweenness  */
  /*************************/


  void exact_betweenness(g_t& G, std::vector<int>& seeds, std::vector<double>& wbtws) {
    int V = G.size();
    std::vector<bool> is_seed(V);


    rep (i, seeds.size()) is_seed[seeds[i]] = true;

    wbtws = std::vector<double>(V);

    rep (s, V) {
      std::vector<double> btws(V);

      std::vector<int> dists, nums, nums_with_seeds;
      bfs(G, s, dists, nums);

      restricted_bfs(G, s, dists, nums_with_seeds, is_seed);

      std::vector<int> out_degree(V);
      rep (u, V) {
        rep (i, G[u].size()) {
          int v = G[u][i];
          if (dists[v] == dists[u] + 1) out_degree[u]++;
        }
      }

      int qh = 0, qt = 0;
      std::vector<int> q(V);
      rep (u, V)
        if (out_degree[u] == 0)
          q[qt++] = u;

      while (qh != qt) {
        int u = q[qh++];
        if (u == s) continue;

        rep (i, G[u].size()) {
          int v = G[u][i];
          if (dists[v] == dists[u] - 1) {
            if (--out_degree[v] == 0) q[qt++] = v;
          } else if (dists[v] == dists[u] + 1) {
            double k = 0;
            k += (double)nums_with_seeds[u] / nums[v];
            if (nums_with_seeds[v] != 0) {
              if (!is_seed[v]) {
                k += btws[v] / nums_with_seeds[v] * nums_with_seeds[u];
              }
            }
            btws[u] += k;
          }
        }
      }
      rep (u, V) {
        if (is_seed[u]) continue;
        wbtws[u] += btws[u];
      }
    }
  }

  void approximate_betweenness(g_t& G, int M, std::vector<int>& seeds, std::vector<double>& btws) {
    int V = G.size();
    whg_t H;
    build_betweenness_hypergraph(G, H, M, seeds);

    btws = std::vector<double>(V);
    rep (i, H.size()) {
      whe_t& whe = H[i];
      rep (j, whe.size()) {
        int v = whe[j].first;
        btws[v] += whe[j].second;
      }
    }
  }

  whe_t make_hyperedge(g_t& G, int s, int t, std::vector<bool>& is_seed, std::unordered_set<int>& vertices_in_whe) {
    whe_t whe;
    std::unordered_map<int, double> btws;
    std::unordered_map<int, int> dists, nums, nums_with_seeds;
    bfs(G, s, dists, nums, vertices_in_whe);

    restricted_bfs(G, s, dists, nums_with_seeds, is_seed, vertices_in_whe);

    std::unordered_set<int> added;
    int qh = 0;
    std::vector<int> q;
    q.push_back(t);
    while (qh != (int)q.size()) {
      int u = q[qh++];
      if (u == s) continue;

      rep (i, G[u].size()) {
        int v = G[u][i];
        if (!vertices_in_whe.count(v)) continue;
        if (!dists.count(v)) continue;

        if (dists[v] == dists[u] - 1) {
          if (!added.count(v)) {
            q.push_back(v);
            added.insert(v);
          }
        } else if (dists[v] == dists[u] + 1) {
          double k = 0;
          k += (double)nums_with_seeds[u] / nums[v];
          if (nums_with_seeds[v] != 0) {
            if (!is_seed[v]) {
              k += btws[v] / nums_with_seeds[v] * nums_with_seeds[u];
            }
          }
          btws[u] += k;
        }
      }
    }
    for (auto it = added.begin(); it != added.end(); ++it) {
      int u = *it;
      whe.push_back(std::make_pair(u, btws[*it]));
    }
    return whe;
  }


  void adaptive_approximate_betweenness(g_t& G, int M, int k, std::vector<int>& seeds, std::vector<double>& btwss) {
    int V = G.size();
    std::vector<bool> is_seed(V);
    seeds.clear();
    btwss.clear();

    whg_t H;
    std::vector<std::pair<int, int > > pairs;
    build_betweenness_hypergraph(G, H, M, seeds, &pairs);

    // preprocess using H
    std::vector<double> degrees(V);
    std::vector<std::vector<int> > vertex_to_heids(V);
    std::vector<std::unordered_set<int> > vertices_in_whes(H.size());
    rep (i, H.size()) {
      whe_t& whe = H[i];
      rep (j, whe.size()) {
        int v = whe[j].first;
        degrees[v] += whe[j].second;
        vertex_to_heids[v].push_back(i);
        vertices_in_whes[i].insert(v);
      }
      vertices_in_whes[i].insert(pairs[i].first);
      vertices_in_whes[i].insert(pairs[i].second);
    }

    std::priority_queue<std::pair<double, int> > pq;
    rep (u, V) pq.push(std::make_pair(degrees[u], u));

    std::vector<double> current_degrees = degrees;
    std::vector<bool> vertex_done(V);


    while (!pq.empty() && (int)seeds.size() < k) {
      double weight = pq.top().first;
      int u = pq.top().second;
      pq.pop();
      if (vertex_done[u]) continue;
      if (fabs(weight - current_degrees[u]) > eps) continue;
      vertex_done[u] = true;

      seeds.push_back(u);
      is_seed[u] = true;
      btwss.push_back(weight);

      rep (i, vertex_to_heids[u].size()) {
        int heid = vertex_to_heids[u][i];
        whe_t& whe = H[heid];

        whe_t new_whe = make_hyperedge(G, pairs[heid].first, pairs[heid].second, is_seed, vertices_in_whes[heid]);
        std::map<int, int> make_sure;
        rep (j, new_whe.size()) make_sure[new_whe[j].first] = j;

        rep (j, whe.size()) {
          int v = whe[j].first;
          if (vertex_done[v]) continue;
          if (whe[j].second != new_whe[make_sure[v]].second) {
            current_degrees[v] -= whe[j].second;
            current_degrees[v] += new_whe[make_sure[v]].second;
            whe[j].second = new_whe[make_sure[v]].second;
          }
          pq.push(std::make_pair(current_degrees[v], v));
        }
      }
    }
  }


  /************************/
  /*  Coverage Centrality */
  /************************/

  bool seed_is_on_the_way(g_t& G, int t, std::vector<int> dists, std::vector<bool>& is_seed) {
    int V = G.size();
    int qh = 0, qt = 0;
    std::vector<int> q(V);
    std::vector<bool> added(V);

    added[t] = true;
    q[qt++] = t;
    if (is_seed[t]) return true;
    while(qh != qt) {
      int u = q[qh++];
      rep (i, G[u].size()) {
        int v = G[u][i];
        if (added[v] == false && dists[v] != -1 && dists[v] == dists[u] - 1) {
          added[v] = true;
          q[qt++] = v;
          if (is_seed[v]) return true;
        }
      }
    }
    return false;
  }



  void exact_coverage(g_t& G, std::vector<int>& seeds, std::vector<int>& btws) {
    int V = G.size();
    std::vector<bool> is_seed(V);
    rep (i, seeds.size()) is_seed[seeds[i]] = true;


    btws = std::vector<int>(V);

    rep (s, V) {
      std::vector<int> dists, nums;
      bfs(G, s, dists, nums);

      rep (t, V) {
        if (seed_is_on_the_way(G, t, dists, is_seed)) continue;

        int qh = 0, qt = 0;
        std::vector<int> q(V);
        std::vector<bool> added(V);

        added[t] = true;
        btws[t]++;
        q[qt++] = t;
        while (qh != qt) {
          int u = q[qh++];
          rep (i, G[u].size()) {
            int v = G[u][i];
            if (added[v] == false && dists[v] == dists[u] - 1) {
              added[v] = true;
              btws[v]++;
              q[qt++] = v;
            }
          }
        }
      }
    }
  }

  void approximate_coverage(g_t& G, int M, std::vector<int>& seeds, std::vector<int>& btws) {
    int V = G.size();
    hg_t H;
    build_coverage_hypergraph(G, H, M, seeds);

    btws = std::vector<int>(V);
    rep (i, H.size()) {
      he_t& he = H[i];
      rep (j, he.size()) {
        int v = he[j];
        ++btws[v];
      }
    }
  }


  void adaptive_approximate_coverage(g_t& G, int M, int k, std::vector<int>& seeds, std::vector<int>& btwss) {
    int V = G.size();
    seeds.clear();
    btwss.clear();

    hg_t H;
    build_coverage_hypergraph(G, H, M, seeds);

    std::vector<double> degrees(V);
    std::vector<std::vector<int> > vertex_to_heids(V);
    rep (i, H.size()) {
      he_t& he = H[i];
      rep (j, he.size()) {
        int v = he[j];
        ++degrees[v];
        vertex_to_heids[v].push_back(i);
      }
    }

    std::priority_queue<std::pair<double, int> > pq;
    rep (u, V) pq.push(std::make_pair(degrees[u], u));

    std::vector<double> current_degrees = degrees;
    std::vector<bool> vertex_done(V);
    std::vector<bool> he_done(H.size());

    while (!pq.empty() && (int)seeds.size() < k) {
      double weight = pq.top().first;
      int u = pq.top().second;
      pq.pop();
      if (vertex_done[u]) continue;
      if (fabs(weight - current_degrees[u]) > eps) continue;
      vertex_done[u] = true;

      seeds.push_back(u);
      btwss.push_back(weight);


      rep (i, vertex_to_heids[u].size()) {
        int heid = vertex_to_heids[u][i];
        if (he_done[heid]) continue;
        he_done[heid] = true;
        he_t& he = H[heid];
        rep (j, he.size()) {
          int v = he[j];
          --current_degrees[v];
          pq.push(std::make_pair(current_degrees[v], v));
        }
      }
    }
  }


  /**************************/
  /*  Building Hypergraphs  */
  /**************************/


  void build_betweenness_hypergraph(g_t& G, whg_t& H, int M, std::vector<int>& seeds, std::vector<std::pair<int, int> >* pairs) {
    int V = G.size();
    std::vector<bool> is_seed(V);
    rep (i, seeds.size()) is_seed[seeds[i]] = true;

    rep (i, M) {
      int s = random() % V, t = random() % V;
      if (pairs) pairs->push_back(std::make_pair(s, t));

      std::vector<int> dists, nums, nums_with_seeds;
      bfs(G, s, dists, nums);

      restricted_bfs(G, s, dists, nums_with_seeds, is_seed);

      int qh = 0, qt = 0;
      std::vector<int> q(V);
      std::vector<bool> added(V);
      std::vector<double> btws(V);
      whe_t whe;
      q[qt++] = t;
      added[t] = true;

      while(qh != qt) {
        int u = q[qh++];
        if (u == s) continue;

        rep (i, G[u].size()) {
          int v = G[u][i];
          if (dists[v] == dists[u] - 1) {
            if (added[v] == false) {
              q[qt++] = v;
              added[v] = true;
            }
          } else if (dists[v] == dists[u] + 1) {
            double k = 0;
            k += (double)nums_with_seeds[u] / nums[v];
            if (nums_with_seeds[v] != 0) {
              if (!is_seed[v]) {
                k += btws[v] / nums_with_seeds[v] * nums_with_seeds[u];
              }
            }
            btws[u] += k;
          }
        }
        whe.push_back(std::make_pair(u, btws[u]));
      }
      H.push_back(whe);
    }
  }


  void build_coverage_hypergraph(g_t& G, hg_t& H, int M, std::vector<int>& seeds) {
    int V = G.size();
    std::vector<bool> is_seed(V);
    rep (i, seeds.size()) is_seed[seeds[i]] = true;

    rep (i, M) {
      int s = random() % V, t = random() % V;

      std::vector<int> dists, nums;
      bfs(G, s, dists, nums);
      if (dists[t] == -1) continue;
      if (seed_is_on_the_way(G, t, dists, is_seed)) continue;

      int qh = 0, qt = 0;
      std::vector<int> q(V);
      std::vector<bool> added(V);
      he_t he;
      q[qt++] = t;
      added[t] = true;
      he.push_back(t);

      while(qh != qt) {
        int u = q[qh++];
        rep (i, G[u].size()) {
          int v = G[u][i];
          if (added[v] == false && dists[v] != -1 && dists[v] == dists[u] - 1) {
            q[qt++] = v;
            added[v] = true;
            he.push_back(v);
          }
        }
      }
      H.push_back(he);
    }
  }


  /*********************/
  /*  Other Functions  */
  /*********************/

  void load_graph(const std::string& fn, g_t& G) {
    std::ifstream ifs(fn.c_str());
    for (std::string line; getline(ifs, line); ) {
      if (line.size() > 0 && line[0] == '#') continue;
      std::istringstream iss(line);
      int u = 0, v = 0;
      iss >> u >> v;
      if (u >= (int)G.size() || v >= (int)G.size()) G.resize(std::max(u, v) + 1);
      G[u].push_back(v);
      G[v].push_back(u);
    }
    rep (u, G.size()) {
      sort(G[u].begin(), G[u].end());
      G[u].erase(unique(G[u].begin(), G[u].end()), G[u].end());
    }
  }



}

