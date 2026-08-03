# Porting notes

## Scope

All nine exported computational routines in `NAMESPACE` are represented. The
R S4 containers `statp` and `adjp` are replaced by Fortran derived types.
Rcpp, Eigen, MASS, and R's `stats` package are not required by the compiled
library.

## Simulation recursions

The C++ recursions in `src/svj.cpp` were translated directly into `lp_path`,
`pvc_path`, and `pv2_path`. High-level simulators generate the same categories
of shocks using a deterministic Park-Miller generator, Box-Muller normals,
Marsaglia-Tsang gamma/chi-square draws, and Poisson rejection sampling.

Random streams are not bit-for-bit identical to R. Identical user seeds are
reproducible within this Fortran implementation.

## Broken upstream jump sampling

Both upstream jump wrappers contain calls of the form:

```text
sapply(x, rnorm, n=1, mean=0)
```

The element supplied by `sapply` conflicts with the named `n` argument to
`rnorm`; therefore these calls do not unambiguously specify a distribution.
The translation implements the apparent intended stochastic models:

- `svj`: if `N` events occur, the jump standard deviation is
  `sqrt(N*sigma1)`.
- `sv1fj`: if `N` events occur, the jump standard deviation is `N`, matching
  the most direct interpretation of passing the count as the missing `sd`
  argument. Multiple events within one very small interval are rare under the
  package defaults.

The actual counts, scales, and realized increments are all returned so callers
can inspect or replace this convention through the low-level recursions.

## Return lengths

The upstream `SV` returns the initial value and therefore has length `n+1`.
`SVJ`, `SV1F`, `SV1FJ`, and `SV2F` remove the initial value and have length
`n`. The port preserves this behavior.

## Probability calculations

The translation includes self-contained implementations of:

- normal CDF and inverse normal CDF;
- regularized incomplete gamma and chi-square CDF;
- regularized incomplete beta and beta CDF;
- R type-7 sample quantiles;
- Pearson correlation and Benjamini-Hochberg adjustment.

For correlation-adjusted pooling, a constant transformed p-value column is
assigned zero off-diagonal correlation instead of propagating R `NA` values.
This gives a deterministic, usable result for degenerate inputs.

## Input validation

The R package relies mainly on downstream errors. The Fortran API checks
matrix shapes, finite values, positive simulation dimensions, valid
correlations, positive CIR parameters, and nondegenerate jump-test samples.
Failures are returned through explicit status values.
