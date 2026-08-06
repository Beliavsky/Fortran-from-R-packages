# sparseIndexTracking modern Fortran

A self-contained modern Fortran translation of the computational code in the
R package `sparseIndexTracking` 0.1.1.

The library computes long-only sparse portfolios that track a benchmark while
satisfying a per-asset upper bound and a full-investment constraint. It uses a
majorization-minimization algorithm with the smooth log-sum approximation to
portfolio cardinality and squared or Huber-type tracking losses.

## Implemented functionality

- Empirical tracking error (`ete`)
- Downside risk (`dr`)
- Huber empirical tracking error (`hete`)
- Huber downside risk (`hdr`)
- Continuation over the upstream sequence of log-sum smoothing parameters
- Two-step squared iterative acceleration with objective backtracking
- KKT projection onto the capped simplex
- Optional initial portfolio and post-fit weight threshold
- Typed fit results with objective, iteration count, cardinality, status, and message
- `spIndexTrack` compatibility entry point and an idiomatic typed API

The R package contains no plotting routine in its computational source. R
objects, `xts` handling, column names, data frames, and package hooks are not
part of the Fortran library.

## Build with FPM

```text
fpm test
fpm run --example example_sparse_index_tracking
```

The package name in `fpm.toml` is `sparse_index_tracking`.

## Build with GNU Make

Checked build with runtime checks:

```text
make check
```

Optimized build:

```text
make release
```

Windows command scripts using `gfortran` and `mingw32-make` are supplied in
`scripts/`.

## Minimal example

```fortran
program fit_tracker
   use sparse_index_tracking, only : dp, sparse_index_fit, &
                                     fit_sparse_index_tracking, sit_success
   implicit none

   real(dp) :: x(100, 5), benchmark(100)
   type(sparse_index_fit) :: fit

   ! Fill x and benchmark with net returns.
   call fit_sparse_index_tracking(x, benchmark, 1.0e-7_dp, fit, &
                                  upper_bound=0.5_dp, measure='ete')
   if (fit%info /= sit_success) error stop trim(fit%message)
   print *, fit%weights
end program fit_tracker
```

A complete deterministic example is in
`example/example_sparse_index_tracking.f90`.

## Main API

Use the facade module:

```fortran
use sparse_index_tracking
```

Important procedures and types are:

- `type(sparse_index_fit)`
- `fit_sparse_index_tracking`
- `sp_index_track`
- `spIndexTrack`
- `tracking_objective`
- `project_capped_simplex`
- `bisection`

The measure names are `ete`, `dr`, `hete`, and `hdr`. The Huber parameter is
required for the last two.

## Numerical and source-compatibility notes

The original algorithm requires `lambda > 0`, although its R documentation
says nonnegative. The Fortran API rejects zero explicitly instead of allowing
a division by zero.

The upstream HDR objective uses R's scalar `&&` operator on vector conditions.
That makes its backtracking objective depend on the sign of the first residual.
The Fortran implementation uses the intended elementwise one-sided Huber loss
by default. Set `source_compatible_hdr_objective=.true.` to reproduce the
upstream objective behavior.

The R code obtains largest symmetric eigenvalues from `eigen()`. This port uses
a deterministic power iteration for positive-semidefinite cross-product
matrices, avoiding a BLAS/LAPACK dependency. Results should be numerically close
but need not be bit-for-bit identical.

The source thresholds small weights and then renormalizes. This can make a
remaining weight exceed `upper_bound`; the Fortran translation preserves that
post-processing behavior.

The source has one global 1000-iteration objective array but can index beyond it
if the limit is reached before the last continuation stage. The Fortran port
stops safely at the requested global iteration limit and returns status
`sit_iteration_limit` with the current feasible portfolio.

## Validation

Five deterministic tests cover the capped-simplex KKT solution, sparse
coefficient recovery, all four tracking measures, robust behavior under an
outlier, exact objective formulas, source-compatible HDR behavior, thresholding,
and invalid inputs.

See `docs/VALIDATION.md` for details.

## License

The upstream package declares `GPL-3`. This translation is distributed under
`GPL-3.0-only`. The complete license is in `LICENSE` and
`license/GPL-3.0.txt`. The unmodified upstream source tree and original archive
are retained in `upstream/` for license and provenance purposes.
