# gRbase

Modern free-form Fortran translation of the portable computational core of the
R package **gRbase 2.0.3**. The translation is intended to live as the
`gRbase/` top-level directory of `Fortran-from-R-packages` and to build with
Fortran Package Manager (FPM).

This is a numerical/algorithmic API, not an emulation of R objects. Graphs use
1-based integer adjacency matrices and tables use explicit integer variable
labels, dimensions, and flattened Fortran/R column-major values.

## Implemented computational areas

- multidimensional cell/entry indexing, slicing, permutations, combinations;
- integer set operations, subset generation, maximal/minimal sets and pairs;
- table construction, permutation, expansion, alignment, marginalization,
  arithmetic, slicing, normalization, and deterministic sampling from supplied
  uniforms;
- row sums, column sums, recycled columnwise products and matrix nonzero indices;
- graph construction and validation, topological sorting, DAG tests,
  moralization, maximum-cardinality search and chordality tests;
- elimination and minimum-cardinality-weight triangulation, recursive thinning
  to an inclusion-minimal triangulation;
- graph queries, separation tests, simplicial nodes, RIP/junction-tree ordering,
  and maximal-prime-subgraph decomposition;
- maximal cliques and connected components through the sibling `igraph`
  translation;
- symmetric-positive-definite inversion and covariance/concentration to partial
  correlation conversion through `rfortran-linalg`.

See `API_COVERAGE.md` for the detailed parity boundary.

## Dependencies

Extract/build this directory beside these existing repository packages:

- `../rfortran-core`
- `../rfortran-linalg`
- `../igraph`

No dependency source is copied into this package. `rfortran-linalg` owns the
pinned `fortran-lapack` dependency, so this package does not link system BLAS or
LAPACK libraries.

## Build and test

```text
fpm build
fpm test
fpm run --example grbase_demo
```

The maintained sources use one real kind, `dp`, imported from `r_kinds` in
`rfortran-core`.

## Small example

```fortran
program demo
  use grbase
  implicit none
  integer, allocatable :: cycle(:, :)
  integer, allocatable :: triangulated(:, :)
  real(dp) :: levels(4)

  cycle = adjacency_from_edges(4, reshape([1, 2, 3, 4, 2, 3, 4, 1], [2, 4]))
  levels = 2.0_dp
  triangulated = minimal_triangulation(cycle, levels)
  print '(a,i0)', 'fill edges: ', (sum(triangulated) - sum(cycle)) / 2
end program demo
```

## Licensing and provenance

The upstream gRbase package is GPL-2-or-later. `COPYING`, `LICENSE`,
`NOTICE.md`, `upstream/DESCRIPTION`, and `upstream/CITATION` preserve the
applicable license and provenance information. The combination-generation
algorithm has an additional historical attribution retained in `NOTICE.md`.
