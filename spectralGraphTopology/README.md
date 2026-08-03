# spectralGraphTopology-fortran

A self-contained modern Fortran translation of the computational code in the R
package **spectralGraphTopology 0.2.3**. The project uses the Fortran Package
Manager (FPM), preserves the original GPL-3 license, and omits only R-specific
packaging, progress-display, and plotting examples.

## Included algorithms

- Graph linear operators: `L`, `A`, `D` and their adjoints/inverses.
- Matrix utilities, structural recovery metrics, pairwise distances, and block
  diagonal construction.
- Spectral graph learning for k-component, cospectral, bipartite, and jointly
  constrained bipartite k-component graphs.
- Constrained Laplacian Rank clustering.
- Smooth graph learning using the Kalofolias primal-dual method.
- Smooth signal-representation graph learning based on Dong et al.
- GLE majorization-minimization and ADMM estimators.
- Combinatorial Graph Laplacian estimation.

The implementation does not require LAPACK, BLAS, R, Eigen, Armadillo, CVXR, or
quadprog. Dense linear algebra, symmetric eigendecomposition, pseudoinversion,
isotonic regression, simplex projection, and nonnegative quadratic optimization
are implemented in Fortran.

## Build

```text
fpm build
fpm test
fpm run
fpm run --example operators_example
fpm run --example spectral_learning_example
fpm run --example smooth_graph_example
fpm run --example laplacian_estimators_example
```

The package version is the valid semantic version `0.2.3`.

## Basic use

```fortran
use spectral_graph_topology, only : dp, graph_result, learn_k_component_graph

type(graph_result) :: graph
real(dp) :: covariance(4,4)

! Fill covariance, then estimate the graph.
call learn_k_component_graph(covariance, graph, k=1, beta=100.0_dp)

if (graph%convergence) then
   print *, graph%laplacian
end if
```

Learning routines are subroutines returning a `graph_result`. Optional arguments
mirror the R defaults where practical. The `use_qp` logical selects the R
package's `"qp"` initialization; the default is its `"naive"` initialization.

See [API.md](API.md), [PORTING.md](PORTING.md), and
[TRANSLATION_COVERAGE.md](TRANSLATION_COVERAGE.md) for details.

## License

The original package declares GPL-3. This translation is distributed under
**GPL-3.0-only**. Selected original metadata and computational source files are retained under
`original/`.
