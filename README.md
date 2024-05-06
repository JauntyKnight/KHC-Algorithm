# Complete Graph Isomorphism Test for Planar Graphs

The repository contains an implementation of a complete graph isomorphism test for planar graphs, 
as per the paper *"Algorithm and Experiments in Testing Planar Graphs for Isomorphism"* (DOI:10.7155/jgaa.00094) by Kukluk, Holder, and Cook, referred to as KHC throughout this repository.


The function `find_planar_code(G)` computes a code for a planar graph `G` that is invariant under isomorphism. Thus, two planar graphs `G1` and `G2` are isomorphic if and only if `find_planar_code(G1) == find_planar_code(G2)`. The algorithm runs in `O(V^2)` time, where `V` is the number of vertices in the graph.


### Key aspects of the algorithm

The KHC algorithm is based on Weinberg's algorithm for building a canonical code for planar triconnected graphs. A quick overview of the KHC algorithm is as follows:

1. Split the graph into biconnected components.
2. Construct a SPQR-tree for each biconnected component, which essentially splits the component into triconnected components.
3. Compute the canonical code of the SPQR-tree, and thus of the biconnected component.
4. Combine the canonical codes of the biconnected components to get the canonical code of the entire graph.

### Key aspects of the implementation

I have created this repository, because to my knowledge, there is no existing open-source implementation of the KHC algorithm (in any programming language).

I have tried to keep the implementation as close to the original paper as possible, with only a few minor optimizations. The code is written in Python, and it uses the `sage` library for graph operations, most importantly for building the SPQR-tree.

### Dependencies

See the file `environment.yml` for the list of dependencies. You can install them using `conda` by running the following command:

```bash
conda env create -f environment.yml
```
