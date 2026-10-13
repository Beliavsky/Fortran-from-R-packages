# HDInterval: modern Fortran numerical translation

Highest-density intervals for draws, density grids and quantile functions, and numerical CDF inversion, translated from **HDInterval 0.2.4**. All maintained Fortran is free form, uses one `dp` kind from `r_kinds`, and all numerical entry points are pure.

## Build and run

Place this directory at `Fortran-from-R-packages/HDInterval`, next to the existing `rfortran-core` directory. With FPM and gfortran, including on Windows:

```text
fpm build
fpm test
fpm run --example --all
```

No separately installed BLAS, LAPACK, R, Python or other runtime is required. The only dependency is `rfortran-core = { path = "../rfortran-core" }`, providing the shared kind and NaN-aware sorting. Dependency code is not copied into this package. The existing `stats` package was checked but could not be built unchanged because its `r_stats_multivariate` module requires the absent `r_matrix_statistics` module. Focused golden-section minimization and CDF bisection are therefore implemented here; see [NOTICE.md](NOTICE.md).

IEEE NaN/infinity behavior is part of the API. Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent assumptions.

## Public API

Import `hdinterval`. All arrays and callbacks use its re-exported `dp` kind. Endpoints are ordered **lower, upper**. Optional status outputs should be checked before using results.

| Procedure | Operation |
|---|---|
| `hdi_samples(samples, interval, cred_mass, status)` | Shortest interval of sorted draws; output `interval(2)`. |
| `hdi_columns(samples, intervals, cred_mass, statuses)` | Independent columns of a draw-by-variable matrix; output shape `(2,n_variables)`. |
| `hdi_density(x, y, intervals, height, cred_mass, allow_split, status, is_hdi, achieved_mass)` | Density-threshold interval(s); allocatable output shape `(2,n_intervals)`. |
| `hdi_function(icdf, interval, cred_mass, tol, status, lower_tail, max_iterations)` | Minimum-width interval from a pure quantile callback. |
| `inverse_cdf(p, cdf, quantile, tol, lower, upper, status, max_iterations, max_extensions)` | Generic scalar CDF inversion. |
| `inverse_cdf(p, cdf, quantiles, tol, lower, upper, status, max_iterations, max_extensions, element_status)` | Vector overload; outputs have `size(p)` elements. |

The specific procedures `inverse_cdf_scalar` and `inverse_cdf_vector` are also public. For callbacks, declare an associated `procedure(scalar_function_interface), pointer`; targets must be pure scalar functions with a single `real(dp), intent(in)` argument and a `real(dp)` result. This pointer interface gives every dummy an explicit `INTENT`, including callback dummies. See [example/interval_demo.f90](example/interval_demo.f90) for working declarations. Distribution parameters can be supplied through a pure wrapper and read-only host association.

### Samples and columns

- `cred_mass` defaults to 0.95 and must be finite and strictly between 0 and 1.
- NaNs, including R-style missing numeric values, are omitted. Infinities are retained.
- The algorithm compares pairs separated by `floor(n * cred_mass)` positions. The chosen interval generally includes `floor(n * cred_mass) + 1` observations. This is upstream's convention, not a quantile interpolation.
- Equal-width ties select the first pair in sorted order. NaN widths caused by equal infinite endpoints are skipped; if all widths are undefined, the sorted range is returned, as in R.
- A singleton has identical endpoints. Empty or all-NaN columns produce NaNs and `hdi_no_interval`; other columns remain available. Input arrays are never changed.
- Chains and posterior draw arrays can be reshaped/combined by the caller before applying these operations. Heterogeneous lists, data frames and MCMC/JAGS objects are not accepted directly.

### Density grids

Coordinates must be finite, increasing and equally spaced, with at least two points. Ordinates must match the coordinate count, be finite and nonnegative, and have a positive maximum. Density normalization need not equal one. A scale-normalized sum prevents overflow for very large density ordinates.

As in upstream, mass is the **sum of grid ordinates**, not a trapezoidal integral. This is appropriate for a regular grid; nonuniform grids are rejected. Spacing is checked with a 128-machine-epsilon coordinate-scaled tolerance. The cutoff retains all ordinates at or above its height, including ties, so attained mass can exceed the requested mass. There is no endpoint interpolation. `height` uses the input's density units.

`allow_split=.true.` returns every selected connected run, including a one-node interval. By default, if several runs are selected, a contiguous credible-interval fallback is returned. It finds the minimum number of grid steps using upstream's strict cumulative-mass comparisons, then averages the tied pairs of endpoint coordinates. `is_hdi` is false for this fallback; the R warning is represented by this flag. When present, `achieved_mass` is the selected ordinate fraction for an HDI, or the **mean cumulative-mass difference of the tied candidate windows** for a fallback. The latter is not a new integration estimate for the averaged endpoints.

A very coarse grid can leave the upstream fallback with no candidate window. R then emits NaN endpoints; this API returns shape `(2,0)` and `hdi_no_interval`. Invalid grids return the same empty shape and `hdi_invalid_data`. Normalizing ordinates changes floating-point summation, so an exactly marginal cumulative cutoff can differ at rounding level from R.

### Quantile callbacks

`hdi_function` searches lower-tail probabilities in `[0,1-cred_mass]`. Its default absolute tail-probability tolerance is `1e-8`, and its iteration limit is 200. A finite median, finite probed quantiles and nonnegative interval widths are required. Endpoints of unbounded distributions need not be finite because the search probes the interior.

The implementation uses **golden-section search**, whereas R uses `stats::optimize` (Brent). Both are bounded unimodal searches, not global optimizers for arbitrary multimodal or irregular width functions. Callbacks must be nondecreasing quantile functions; local checks cannot prove that property everywhere. For nonunique minima, such as uniform distributions, interval locations may differ from R while widths agree. Boundaries are approached to the specified tolerance, not snapped exactly to zero or one. An exhausted search returns NaNs and `hdi_no_convergence`.

### Numerical CDF inversion

Probabilities must be finite and in `[0,1]`. **All values are clipped to `[1e-6,1-1e-6]`**, exactly as in upstream, so probabilities zero and one yield finite interior quantiles rather than mathematical support endpoints. `log.p` and `lower.tail` are not implemented upstream and have no Fortran options.

The default initial bracket is `[-10,10]`. Missing brackets are extended symmetrically with doubling steps, up to 100 extensions by default, followed by bisection with up to 1000 refinements. The default absolute coordinate tolerance is `epsilon(1.0_dp)**0.25`, about `1.22e-4`. A root can also terminate at an exact probability match or when there is no representable interior coordinate. This differs from R's Brent-style `uniroot` iteration and extension path; it does not claim identical final bits.

A callback must be a continuous, nondecreasing CDF valid over the probed domain. Every sampled value must be finite and in `[0,1]`; the bracket endpoints and bisection probes are checked for monotonicity. These checks cannot certify an arbitrary callback globally. For a discontinuous CDF, bisection may approach a jump without attaining the requested probability; generalized discrete quantiles are not this routine's contract. A constant callback can match one probability without describing a proper distribution. Callers remain responsible for supplying a valid CDF.

Callbacks defined only on a bounded domain should be supplied explicit `lower` and `upper` bounds that already bracket the requested probability, normally with `max_extensions=0`. An unavailable/invalid callback value is reported as `hdi_invalid_callback`, never silently clamped. Scalar failures produce NaN. Vector inversions retain successful elements, provide optional per-element statuses, and return the first failure as the aggregate status. An empty vector succeeds without calling or validating a callback.

## Status constants

| Value | Name | Meaning |
|---:|---|---|
| 0 | `hdi_ok` | Computation succeeded. |
| 1 | `hdi_invalid_mass` | Credible mass is not finite or not strictly between 0 and 1. |
| 2 | `hdi_invalid_data` | Invalid probability, shape, density ordinate or grid. |
| 3 | `hdi_invalid_callback` | Nonfinite, out-of-range or observed nonmonotone callback result. |
| 4 | `hdi_no_interval` | No usable sample or no contiguous density fallback window. |
| 5 | `hdi_no_convergence` | Bracket extension or optimization/root refinement limit exhausted. |
| 6 | `hdi_invalid_options` | Invalid tolerance, bracket or iteration limit. |
| 7 | `hdi_unassociated_callback` | Null procedure pointer. |

## Translation coverage

**Status: substantial. 2 of 2 (100%)** exported computational R functions are mapped. The fraction measures mapped computational functions, not complete R compatibility. `NAMESPACE` exports only `hdi` and `inverseCDF`; the S3 methods are registered but not separately exported, and are accounted for under `hdi` rather than inflating the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `hdi` | `hdinterval` | `hdi_samples`, `hdi_columns`, `hdi_density`, `hdi_function` | substantial |
| `inverseCDF` | `hdinterval` | `inverse_cdf`, `inverse_cdf_scalar`, `inverse_cdf_vector` | substantial |

No exported computational function is left unmapped. R dispatch, class attributes, result names, package-specific object extraction, warnings, and `...` are omitted. [API_COVERAGE.md](API_COVERAGE.md) inventories all registered methods and numerical helpers. The manifest contains matching machine-readable coverage.

## Verification and provenance

`fpm test` runs deterministic tests against actual HDInterval 0.2.4 R outputs and independent analytic distribution identities. Cases include sample rounding/ties, NaNs and infinities, column-wise missingness, density cutoff ties and scale invariance, multimodal splitting and fallback, normal/logistic/exponential/Beta/Cauchy quantile intervals, clipped inverse probabilities, far-away roots, malformed callbacks and explicit iteration failures.

The retained [reference/generate_reference.R](reference/generate_reference.R) sources the unmodified upstream R files to regenerate [reference/expected.txt](reference/expected.txt) and [test/r_reference_values.f90](test/r_reference_values.f90). Run it from the package directory with `Rscript reference/generate_reference.R`. R is needed only for regeneration. The fixtures were generated with R 4.5.0. Tests use exact comparisons for discrete sample endpoints, approximately `1e-14` for grid results, `1e-7` to `3e-7` for quantile-interval comparisons across different optimizers, and `2e-4` for default-tolerance root comparisons. Tight-tolerance logistic inversion is independently checked within `2e-9`.

Licensed under **GPL-3.0-only**, preserving the upstream GPL-3 designation. See [LICENSE](LICENSE), [NOTICE.md](NOTICE.md), and [PROVENANCE.md](PROVENANCE.md). The sibling dependency retains its own MIT license.
