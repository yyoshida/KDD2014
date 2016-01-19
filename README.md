# Adaptive Betweenness/Coverage Centrality
A library for computing adaptive betweenness/coverage centrality.

## Usage

From CUI Interface

    $ ./waf configure
    $ ./waf
    $ ./build/default/test -G=sample-graph.txt
*  Execute waf to build a test program.
*  Execute test to compute vertex centralities of the sample graph.

See the source code to figure out how to use the library from your program.
In the graph file, each line should contain two vertices (see sample-graph.txt).
Vertices should be numbered from zero.

## Reference
* Yuichi Yoshida. 2014. Almost linear-time algorithms for adaptive betweenness centrality using hypergraph sketches. In Proceedings of the 20th ACM SIGKDD international conference on Knowledge discovery and data mining (KDD'14). ACM, New York, NY, USA, 1416-1425. DOI=http://dx.doi.org/10.1145/2623330.2623626



