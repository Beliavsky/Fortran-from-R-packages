# API coverage

This document summarizes computational parity with the attached gRbase 2.0.3
snapshot. Names are Fortran-native rather than R/Rcpp ABI names.

## Arrays, cells and combinations

| Upstream area | Fortran coverage |
| --- | --- |
| `cell2entry_`, `entry2cell_`, `make_plevels_` | `cell_to_entry`, `entry_to_cell`, `make_plevels` |
| `next_cell_`, slice cell iteration | `next_cell`, `next_cell_slice`, `slice_to_entries` |
| permuted cell entries | `cell_to_entry_perm`, `perm_cell_entries` |
| `fastcombn`, Bristol combination primitive | `choose_integer`, `combinations`, `next_combination_mask` |
| row/column reductions | `row_sums`, `column_sums` |
| `colwiseProd` | `columnwise_product` with cyclic weight recycling |
| `which_matrix_index__` | `matrix_nonzero_indices` with 1-based `(row,column)` output |

## Set operations

Implemented: sorting/deduplication, union, intersection, difference, subset
checks, containment, all subsets, maximal sets, minimal sets and all unordered
pairs. Integer labels replace R character-vector dispatch. The compiled
`filter_maximal_vectors_` behavior is covered by `maximal_sets`.

## Multidimensional tables

Implemented using `table_t` (`var`, `dim`, flattened `value`):

- construction and validation;
- permutation and alignment;
- expansion with replicate/divide/zero-fill modes;
- marginalization and slicing;
- add/subtract/multiply/divide, zero-safe division, and list add/multiply folds;
- equality independent of variable order;
- global and first-dimension conditional normalization;
- deterministic weighted sampling through `table_sample_from_uniforms`.

The sampler accepts caller-supplied uniforms rather than owning an RNG. This
keeps the computational sampling transform while excluding R's RNG state,
`set.seed`, and R `sample` interface.

## Graph algorithms

Implemented directly on integer adjacency matrices:

- graph construction and adjacency validation;
- directed topological sort and DAG test;
- moralization;
- maximum-cardinality search and chordality test;
- elimination-order triangulation;
- minimum-cardinality-weight heuristic triangulation;
- recursive thinning to an inclusion-minimal triangulation;
- neighbors, parents, children, ancestors and descendants;
- complete-set and simplicial-node queries;
- set separation;
- RIP/junction-tree ordering;
- maximal-prime-subgraph decomposition.

`maximal_cliques_adjacency` and `connected_components_adjacency` deliberately
reuse the sibling `igraph` translation rather than duplicate general graph
algorithms.

## Linear algebra/statistics

Implemented:

- SPD matrix inversion (`inverse_spd`);
- concentration-to-partial-correlation conversion;
- covariance-to-partial-correlation conversion.

These call the sibling `rfortran-linalg` API. No system BLAS/LAPACK link is
introduced by gRbase.

## Intentionally excluded / not parity targets

- R/Rcpp registration and `.Call` wrappers;
- S3/S4 methods, formulas, dimnames and R class coercions;
- dense/sparse `Matrix` class conversion machinery;
- igraph object conversion wrappers already represented by the sibling Fortran
  graph package;
- plotting and interactive helpers;
- datasets, documentation examples tied to R objects, and package-management
  utilities;
- hierarchical log-linear-model convenience interfaces that mostly orchestrate
  R formula/model objects rather than expose an independent numerical kernel;
- internal character/string specializations where integer variable labels give
  the native Fortran equivalent.

## Fidelity notes

Table values use the same first-index-fastest column-major convention as R and
Fortran. Graph vertices and table levels are 1-based. Expansion order follows
upstream gRbase: variables unique to the first table are fastest, followed by
variables from the auxiliary/second table.
