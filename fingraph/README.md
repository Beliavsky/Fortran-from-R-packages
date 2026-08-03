# fingraph-fortran

A self-contained modern Fortran translation of the computational code in the R
package **fingraph 0.1.0**. The project uses the Fortran Package Manager (FPM)
and implements the package's ADMM estimators for financial graph learning under
Gaussian and Student-t assumptions.

The graph operators, pseudoinverse, eigensolver, and constrained initialization
were adapted from the earlier `spectralGraphTopology-fortran` translation. The
portable random-number and Student-t simulation support used by the tests and
examples was adapted from `fitHeavyTail-fortran`.

## Included algorithms

- `learn_connected_graph`: connected combinatorial Laplacian estimation from a
  covariance matrix under a Gaussian model.
- `learn_regular_heavytail_graph`: connected graph estimation from observations
  under Gaussian or multivariate Student-t assumptions.
- `learn_kcomp_heavytail_graph`: graph estimation with a prescribed number of
  connected components under Gaussian or Student-t assumptions.
- Student-t observation weights and both augmented-Lagrangian evaluators.
- Laplacian/adjacency operators, adjoints, pseudoinverse, Jacobi symmetric
  eigendecomposition, nonnegative QP initialization, graph metrics, and portable
  Gaussian/Student-t simulation support.

The implementation is dense and dependency-free. It is intended for research,
validation, and moderate graph sizes. Large financial networks will benefit from
future BLAS/LAPACK and sparse-matrix backends.

## Build

```text
fpm build
fpm test
fpm run
fpm run --example connected_covariance_example
fpm run --example student_t_graph_example
fpm run --example k_component_graph_example
```

The FPM version is the valid semantic version `0.1.0`.

A direct GNU Fortran validation script is also included:

```text
./run_gfortran_tests.sh
```

On Windows:

```text
run_gfortran_tests.bat
```

## Basic use

```fortran
use fingraph, only : dp, fingraph_result, learn_regular_heavytail_graph

type(fingraph_result) :: graph
real(dp) :: returns(1000,20)

! Fill returns with centered observations.
call learn_regular_heavytail_graph(returns,graph, &
   heavy_type='student',nu=5.0_dp,rho=1.0_dp)

if (graph%convergence) then
   print *, graph%laplacian
end if
```

R lists are represented by `type(fingraph_result)`. Diagnostic vectors are
allocated to exactly the number of completed iterations. The initial graph can
be supplied explicitly with `initial_weights`, or selected through
`initialization='naive'` or `initialization='qp'`.

For the workflow shown in the original README, estimate `nu` with the separate
`fitHeavyTail-fortran` project and pass the resulting degrees of freedom to the
Student-t Fingraph estimator.

See [API.md](API.md), [PORTING.md](PORTING.md),
[TRANSLATION_COVERAGE.md](TRANSLATION_COVERAGE.md), and
[TESTING.md](TESTING.md).

## License

The R package metadata declares GPL-3, while its repository also contains an
MIT-form license file. This translation incorporates GPL-3.0-only code adapted
from `spectralGraphTopology-fortran` and `fitHeavyTail-fortran`, so the combined
Fortran project is distributed under **GPL-3.0-only**. The original MIT file is
retained verbatim as `original/LICENSE-MIT`; see [NOTICE.md](NOTICE.md).
