# Complete Graph Isomorphism Test for Planar Graphs

The repository contains an implementation of a complete graph isomorphism test for planar graphs, 
as per the paper *"Algorithm and Experiments in Testing Planar Graphs for Isomorphism"* (DOI:10.7155/jgaa.00094) by Kukluk, Holder, and Cook, referred to as KHC throughout this repository.


The function `find_planar_code(G)` computes a code for a planar graph `G` that is invariant under isomorphism. Thus, two planar graphs `G1` and `G2` are isomorphic if and only if `find_planar_code(G1) == find_planar_code(G2)`. The algorithm runs in `O(V^2)` time, where `V` is the number of vertices in the graph.


### Key aspects of the algorithm:

Since Weinberg provided a complete isomorphism test for planar triply connected graphs (DOI: 10.1109/TCT.1966.1082573), the KHC algorithm can be thought of building the code for a graph in a bottom-up fashion:

1. Compute the canonical codes of the triply connected components of the graph.
2. Compute the canonical codes of the biconnected components of the graph, using the canonical codes of the triply connected components.
3. Compute the canonical code of the graph, using the canonical codes of the biconnected components.
