# igraph-fortran 0.3.0

A native modern Fortran/FPM translation of a substantial computational core of
R `igraph` 2.3.3 / its bundled C igraph library.

This is **not** an `iso_c_binding` wrapper. Algorithms in `src/` are native
Fortran 2018. The complete supplied R/C source tree is retained under
`original/igraph-master/` for provenance and future expansion.

## Highlights in v0.3.0

v0.3.0 keeps the complete v0.1/v0.2 API and adds five more algorithm families.

### Canonical labeling and automorphisms

- `canonical_labeling`
- `canonical_code`
- `canonical_form`
- `automorphism_count`
- `automorphisms`
- optional integer vertex colors

The native implementation uses equitable color refinement followed by
individualization/backtracking. `automorphisms` can return explicit maps, while
`automorphism_count` counts the complete group without storing all maps.

This is an exact topology implementation for simple graphs, but it does not
attempt to reproduce the internal bliss data structures used by C igraph.

### Community detection

- `cluster_fast_greedy`
- `cluster_leiden`
- retained weighted `cluster_louvain`
- retained label propagation

`cluster_fast_greedy` performs exact greedy modularity merges and retains the
best partition on the merge path. The Leiden routine uses local modularity
moves, nested within-community refinement, and repeated graph aggregation.
It implements the major Leiden architecture natively, but is not a bit-for-bit
copy of igraph's highly optimized Leiden C implementation.

### Flow and cuts

- `gomory_hu_tree`
- `min_cost_flow`
- retained Dinic maximum flow / minimum cut

`gomory_hu_tree` uses n-1 native maximum-flow solves and returns the weighted
cut tree plus parent/cut arrays. `min_cost_flow` is a successive shortest
augmenting-path residual solver using Bellman-Ford, so negative edge costs are
supported when no reachable negative-cost residual cycle makes the problem
ill-defined.

### Larger motifs / RAND-ESU

- `motif_census` for connected induced motifs of size 2..6
- `randesu_motif_census` with per-level branch-cut probabilities
- canonical motif codes that aggregate isomorphic occurrences
- retained dyad and exact 16-class triad census

The exact enumerator uses the ESU connected-subgraph expansion. Motif
canonicalization is exact by permutation for the supported small motif sizes.

### Products and graph algebra

- Cartesian product
- tensor/direct product
- strong product
- graph union
- graph intersection
- graph difference
- graph join
- retained induced/edge subgraph, transpose, complement, disjoint union and
  line graph

## v0.2 functionality retained

- VF2-style graph/subgraph isomorphism with colors and induced matching
- all/maximal/largest cliques and clique number
- dyad and Davis-Leinhardt 16-class triad census
- exact global edge and vertex connectivity
- source/sink minimum cut and biconnectivity
- weighted modularity and multilevel Louvain
- graph transformations

## v0.1 functionality retained

- edge-list plus outgoing/incoming CSR graph storage
- real edge weights, directed/undirected graphs
- graph construction and simplification
- BFS/DFS/topological sorting
- BFS, Dijkstra, Bellman-Ford and Floyd-Warshall paths
- weak/strong components, articulation points and bridges
- closeness, Brandes betweenness, PageRank and eigenvector centrality
- density, reciprocity, triangles/transitivity and assortativity
- Kruskal MST/forest
- Dinic maximum flow
- Hopcroft-Karp bipartite matching
- Erdos-Renyi, ring, star and tree generators

## Public module

```fortran
use igraph
```

Example:

```fortran
program demo
  use igraph
  implicit none
  type(graph_t) :: g
  type(community_result_t) :: c
  type(canonical_result_t) :: canon

  g = ring_graph(5)
  canon = canonical_labeling(g)
  print *, canon%permutation
  print *, automorphism_count(g) ! 10

  c = cluster_leiden(g, seed=2026_i8)
  print *, c%membership
end program demo
```

## Build

```text
fpm build
fpm test
fpm run --example v3_algorithms
```

The release was additionally validated directly with GNU Fortran using:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

and with `-O2`.

## Current scope

R igraph exposes a very large API, so v0.3.0 remains an expanding native
translation rather than a claim that every specialized igraph routine is
already present. Important deferred areas include:

- optimized bliss internals and large-graph canonical labeling
- Infomap, walktrap, spinglass and edge-betweenness communities
- full graphlet/orbit and larger motif APIs
- k-shortest paths and exhaustive simple-path families
- weighted Brandes betweenness and Johnson all-pairs shortest paths
- spectral sparse eigensolvers / ARPACK methods
- advanced random graph models and degree-sequence samplers
- GraphML/GML/Pajek/DL readers and writers
- layout algorithms
- R attributes, S3/vctrs objects, plotting and callback infrastructure

## Licensing

The R package and bundled igraph core are GPL version 2 or later. This
translation is distributed under GPL-2.0-or-later. See `COPYING` and
`LICENSES.md`.
