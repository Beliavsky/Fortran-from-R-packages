# Validation

Validated with GNU Fortran using two configurations.

Strict debug:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

Optimized:

```text
-std=f2018 -O2
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

Fourteen regression executables pass in both configurations:

1. `test_graph_core`
2. `test_paths_components`
3. `test_centrality`
4. `test_structural`
5. `test_optimization`
6. `test_generators_community`
7. `test_vf2_cliques`
8. `test_connectivity_transform`
9. `test_louvain`
10. `test_motifs`
11. `test_canonical_v3`
12. `test_flow_v3`
13. `test_community_v3`
14. `test_motifs_products_v3`

New v0.3 checks include:

- C4 automorphism group size 8 and K4 group size 24
- colored C4 automorphism group size 4
- canonical-code invariance under a nontrivial vertex relabeling
- explicit automorphism-map enumeration
- min-cost flow sends demand 3 at minimum cost 7 on the regression network
- every Gomory-Hu tree path minimum equals a direct max-flow value for all
  vertex pairs in the test graph
- weighted fast-greedy and Leiden recover the two-triangle partition with
  modularity 0.48360655737704927
- exact ESU finds all five K4 occurrences in K5
- RAND-ESU with zero cuts agrees exactly with ESU; probability-one cuts return
  no occurrences
- Cartesian C3 x C4 has 12 vertices / 24 edges
- tensor C3 x C4 has 12 vertices / 24 edges
- strong C3 x C4 has 12 vertices / 48 edges
- graph join and set-operation edge counts

Representative v0.1/v0.2 checks remain:

- weighted shortest-path distance = 6
- MST total weight = 6
- maximum flow = 5
- bipartite matching size = 3
- cycle C5 edge/vertex connectivity = 2
- K4 edge/vertex connectivity = 3
- VF2 relabeling and subgraph matching
- all 16 Davis-Leinhardt triad classes
- weighted Louvain two-community modularity = 0.48360655737704927

All four examples compile and run under the optimized configuration.
All Fortran source lines are at most 132 columns.
