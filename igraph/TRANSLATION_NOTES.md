# Translation notes

## Source

- R package: `igraph` 2.3.3
- bundled C core: igraph under `src/vendor/cigraph`
- supplied archive: `igraph-master.zip`

The complete supplied tree is retained under `original/igraph-master/`.

## Design

The R package is a high-level interface over a large C library. This project
uses a native Fortran graph representation instead of calling that C library.
`graph_t` stores one canonical edge list and builds outgoing/incoming CSR-like
adjacency arrays. Undirected edges are represented in both adjacency directions
but retain one canonical edge id and weight.

The translation deliberately avoids relying on short-circuit evaluation of
Fortran `.and.` and `.or.` when an array index is protected by a condition.

## v0.3 algorithm mapping

- canonical labeling: equitable color refinement plus recursive
  individualization/backtracking; the lexicographically minimum adjacency code
  defines the canonical order
- automorphisms: exact self-isomorphism enumeration with degree/color/mapped
  adjacency pruning
- fast greedy: repeated weighted modularity-delta community merges while
  retaining the best partition on the merge path
- Leiden-style communities: local moving, within-parent refinement, aggregation,
  then another local-moving level initialized from the coarse partition
- Gomory-Hu: standard parent/cut-tree algorithm with n-1 Dinic max-flow solves
- minimum-cost flow: residual network with successive shortest augmenting paths
  and Bellman-Ford shortest paths
- ESU motifs: connected induced-subgraph expansion with exclusive-neighborhood
  extension; RAND-ESU uses per-level probabilistic branch cuts
- motif classification: exact canonical adjacency bit mask over all motif vertex
  permutations for motif sizes up to six
- graph products: direct adjacency-rule construction for Cartesian, tensor and
  strong products

The v0.3 canonical-labeling implementation is algorithmically exact for simple
unweighted topology but does not reproduce bliss's specialized partition stack,
pruning invariants or multigraph edge-color callbacks. Likewise, `cluster_leiden`
implements the major Leiden local-move/refinement/aggregation architecture but
is not intended to be random-number-for-random-number identical to C igraph.

## v0.2 mapping retained

- VF2-style graph and subgraph matching
- Bron-Kerbosch maximal cliques
- Davis-Leinhardt triad classification
- exact edge/vertex connectivity through max-flow and vertex splitting
- weighted Louvain aggregation
- induced/edge/complement/transpose/line/disjoint-union transformations

## v0.1 mapping retained

- BFS/DFS and topological traversal
- weak components and Kosaraju SCCs
- BFS/Dijkstra/Bellman-Ford/Floyd-Warshall paths
- Brandes unweighted betweenness
- PageRank/eigenvector power iteration
- Tarjan articulation points and bridges
- Kruskal MST
- Dinic maximum flow
- Hopcroft-Karp bipartite matching
- randomized label propagation

## Semantic differences from R igraph

- Vertices are integers `1..n`; R names and arbitrary attributes are outside the
  native numeric graph type.
- Betweenness remains unweighted in v0.3.0.
- Canonical labeling treats adjacency topology; edge weights are carried through
  `canonical_form` but do not participate in canonical partition refinement.
- Motif census v0.3 supports connected induced motifs of sizes 2..6.
- `min_cost_flow` currently targets directed networks.
- Layout, plotting and R object infrastructure are intentionally excluded.
