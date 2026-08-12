# Changelog

## 0.3.0

- Added exact individualization/refinement canonical labeling and canonical forms.
- Added exact automorphism counting and optional explicit automorphism maps.
- Added weighted fast-greedy modularity communities.
- Added a native Leiden-style local-move/refinement/aggregation community solver.
- Added Gomory-Hu cut trees using the native Dinic max-flow backend.
- Added successive-shortest-path minimum-cost flow with residual Bellman-Ford.
- Added connected motif census for sizes 2..6 using ESU.
- Added RAND-ESU branch-cut sampling and canonical motif codes.
- Added Cartesian, tensor and strong graph products.
- Added graph union, intersection, difference and join.
- Added four regression executables and the `v3_algorithms` example.
- Retained every v0.1/v0.2 API and regression test.

## 0.2.0

- Added VF2-style graph and subgraph isomorphism with vertex mappings, optional
  vertex colors and induced matching.
- Added all/maximal/largest clique enumeration and clique number.
- Added Davis-Leinhardt dyad and 16-class triad census.
- Added exact global edge connectivity and vertex connectivity plus source/sink
  minimum-cut and biconnectivity tests.
- Upgraded modularity to use edge weights and added multilevel Louvain
  community detection with graph aggregation.
- Added induced/edge subgraphs, transpose, complement, disjoint union and line
  graph transformations.
- Added four regression executables and the `advanced_algorithms` example.
- Retained all v0.1.0 APIs and tests.

## 0.1.0

- Initial modern Fortran/FPM translation.
- Added native graph storage with edge list and outgoing/incoming CSR arrays.
- Added construction, traversal, shortest paths, components, centralities,
  structural statistics, MST, maximum flow, bipartite matching, generators,
  label propagation and modularity.
- Added six regression executables and two examples.
- Retained the complete original R/C source tree for provenance.
