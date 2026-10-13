# BurStFin for modern Fortran

A free-form Fortran/FPM translation of the computational routines in
**BurStFin 1.3** by Burns Statistics Ltd. It provides statistical-factor and
Ledoit–Wolf equal-correlation covariance estimators, benchmark covariance
transformations, alpha proxies and sliding time weights.

## Build and test

Extract this directory as `BurStFin` at the root of
https://github.com/Beliavsky/Fortran-from-R-packages, beside the existing
`rfortran-core` and `rfortran-linalg` directories. Run:

```sh
fpm build
fpm test
fpm run --example --all
```

The two siblings are dependencies, not included copies. FPM obtains the
unchanged linear-algebra sibling's pinned pure-Fortran LAPACK dependency.
No separately installed BLAS or LAPACK library is required. The manifest and
sources support GNU Fortran/FPM on Windows; validation here uses Linux GNU
Fortran 14.2.0 and FPM 0.13.0, not a native Windows execution.
Do not use `-ffast-math`, `-Ofast`, or finite-only assumptions: missing-value
behavior depends on IEEE NaNs. All maintained reals use `dp` from `r_kinds`.

## Translation coverage

Package status: **substantial**. **6 of 6 (100%)** exported computational R
functions are mapped. The fraction measures mapped computational functions,
not complete R compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| alpha.proxy | burstfin_portfolio | alpha_proxy, alpha_proxy_grid | substantial |
| factor.model.stat | burstfin_factor | factor_model_stat, factor_covariance | substantial |
| slideWeight | burstfin_portfolio | slide_weight | complete |
| var.add.benchmark | burstfin_portfolio | var_add_benchmark | substantial |
| var.relative.benchmark | burstfin_portfolio | var_relative_benchmark | substantial |
| var.shrink.eqcor | burstfin_shrink | var_shrink_eqcor | substantial |

`partial.rainbow` and `threeDarr` are excluded presentation/R-array-interface
helpers. The registered, unexported fitted method is covered numerically by
`factor_covariance`, without adding it to the denominator. See
[API_COVERAGE.md](API_COVERAGE.md) and the exact manifest mappings.

## Public API

Use `use burstfin` for all public procedures, `dp`, `alpha_surface`, and
`statistical_factor_model`, or import their narrower modules.
Arrays use observations in rows and assets in columns. Asset identity is
explicitly positional. No name matching, timestamps, call objects, warning
printing, graphics or R recycling is implemented. All scalar options have
FORD argument documentation in `src/`.

### Statistical factor model

```fortran
call factor_model_stat(returns, model, status)
call factor_covariance(model, covariance, status)
call factor_covariance(model, systematic, status, output='systematic')
call factor_covariance(model, specific, status, output='specific')
```

`model` is a `statistical_factor_model` with allocatable `loadings`,
`uniquenesses`, `sdev`, `cumulative_variance_fraction`, and logical
`constant_assets`. It also records `iterations` and `converged` as useful
extra diagnostics. The default/full covariance is the sum of systematic and
specific covariance. Factor loadings are in standardized units; covariance
is in squared original return units.

All numerical controls are supported: `weights`, logical `center`, explicit
`center_values`, `frac_var`, `iter_max`, integer-vector `nfac_miss`, `full_min`,
`reg_min`, `sd_min`, `quan_sd`, `tol`, `zero_load`, integer `range_factors(2)`,
`constant_returns_okay`, `specific_floor`, and exact string `floor_type`
(`quantile` or `fraction`). Defaults match R. `center_values` overrides
`center`; integer count options replace R's numeric-to-integer coercion.

All nonfinite observations are missing. Entirely missing rows are removed,
and supplied original-length weights are subset accordingly. Weighted
centers deliberately retain the upstream non-reweighting for missing cells.
Complete-asset SVD, principal-factor iterations, missing-asset score
regressions, volatility quantiles, loading rescaling and uniqueness floors
follow the source. In particular, the missing-asset regression intentionally
uses the asset as the predictor of scores, then divides its coefficients by
squared original singular values; it is not silently replaced with a
conventional reverse regression.

### Equal-correlation shrinkage

```fortran
call var_shrink_eqcor(returns, covariance, shrinkage, status)
```

`shrinkage` receives the actual clipped intensity. Controls are `weights`,
explicit `shrink`, logical `center`, explicit `center_values`, `vol_shrink`,
integer `sd_min`, `quan_sd`, `tol`, and `compatible`. Defaults match R.
Zero-weight rows are removed. NaNs are pairwise missing; infinities are
rejected, unlike the factor routine which treats them as missing. Pairwise
sample covariance is computed through the existing core API. The
`compatible` option rescales by `(n-1)/n`, including after row removal.
Volatility shrinkage is applied to the target, while mean correlation uses
the unshrunk sample diagonal. The automatic fourth-moment estimate and the
conditional eigenvalue floor follow the upstream algorithm.

For both covariance estimators, absent weights mean a linear sequence from
0.5 to 1.5. A one-element positive `weights` array means equal weights, the
same numerical result as R's `weights=NULL` or positive scalar. All-zero or
negative weights are invalid.

### Portfolio utilities

- `alpha_proxy(weight,vol_man,vol_bench,vol_other,cor_man,cor_bench)` is pure
  elemental, with all six arguments optional and defaulting to 0.2. It returns
  the upstream proxy in squared-volatility basis-point units (factor 10000).
  Scalar broadcasting and conformable arrays use Fortran elemental rules.
  Invalid domains/nonfinite inputs return a quiet NaN.
- `alpha_proxy_grid(...six real vectors...,surface,status,transpose_output)`
  supports at most two nonscalar inputs. With two, both axes are sorted in
  upstream argument order; `transpose_output=.true.` exchanges them. The
  `alpha_surface` result includes `x`, `y`, `z`, and integer axis-argument
  positions 1 through 6. For a scalar or one-axis result, the unused axis has
  length one and coordinate zero; one varying input preserves its order.
- `slide_weight(n,weights,status,fractions,observations,locations)` supports
  all three coordinate forms. Observations and locations are exclusive;
  either overrides fractions. Endpoints use R's nearest-even rounding,
  then the upstream clipped ramp. The default endpoint at `n` has weight
  `n/(n+1)`, not one.
- `var_add_benchmark(variance,weights,result,status,sum_to_one,asset_indices)`
  accepts a matrix or rank-three cube. The default normalizes the sum of
  absolute weights to one when more than `1e-10` away. Optional unique
  one-based `asset_indices` map sparse/reordered holdings to input assets.
  With normalization disabled, zero holdings are valid.
- `var_relative_benchmark(variance,benchmark,result,status)` accepts a matrix
  or cube and removes the one-based benchmark index. Upstream use of the
  benchmark column on both sides, even for nonsymmetric input, is retained.
  NaNs propagate in this transformation, as in R.

## Error handling and deliberate numerical differences

Status codes are 0 success, 1 invalid input/shape/option, 2 insufficient data,
3 degenerate/invalid numerical state, and 4 shared solver failure. Utilities
use 0/1 only. Do not consume outputs after a nonzero status; allocated pieces
may be present but are not a valid result. Shared eigensolver/SVD routines
are not declared pure, so their callers are ordinary subroutines; all other
maintained computational procedures are pure or pure elemental.

- Covariance inputs are internally scaled before quadratic/fourth moments.
  Explicit shrinkage uses individual asset scales to support very different
  asset magnitudes; automatic shrinkage uses one common scale to preserve
  the fourth-moment estimate.
  This preserves the mathematical algorithm but changes floating-point
  rounding and avoids intermediate underflow/overflow. Tests include global
  rescaling by `1e140` and `1e-140`. Final covariance must still be representable. Automatic shrinkage can reject
  a mixed-scale dataset if global normalization underflows an asset variance;
  an explicit shrink intensity avoids that limitation.
- SVD/eigenvector signs and bases within exactly repeated eigenvalue spaces
  can differ by LAPACK implementation. Compare reconstructed covariances,
  communalities and explained variance rather than arbitrary column signs.
- The R SVD failure retry is not reproduced; shared solver failure returns 4.
  The unchanged shared LAPACK dependency can print a `path 5t` diagnostic
  on some wide-matrix SVD paths; this is a dependency side effect.
- `range_factors=[0,0]` produces a coherent zero-factor model. The R code's
  `1:0` subscript construction can fail; this documented extension avoids it.
- Automatic shrinkage uses intensity zero when the squared sample-target difference is zero
  or below the smallest normal number in the internally scaled units. Then any intensity gives the same matrix. Upstream
  can compute `0/0`; the translated result is defined.
- Singular missing-score regressions, an entirely constant factor universe,
  and zero reliable shrinkage variances return 3 rather than silently making
  undefined NaNs. Permitted individual constant factor assets have zero
  standard deviations and covariance rows/columns.
- Invalid/empty dimensions, infinite controls, and invalid quantile ranges
  are rejected explicitly. R partial string matching, numeric coercion,
  recycling and warning side effects are not reproduced. Fractional
  uniqueness floors above one remain permitted, matching upstream behavior.

## Verification and licenses

The self-contained Fortran test suite includes actual R 4.5.0/BurStFin 1.3
fixtures, analytic rank-one/benchmark/alpha cases, all covariance outputs,
missingness and incomplete-asset regression, time weights and explicit
centers, scale homogeneity, degeneracy, nonfinite values and boundaries.
R is only required to regenerate fixtures with
`Rscript reference/generate_fixtures.R`; it is not needed for FPM builds/tests.
See [reference/README.md](reference/README.md) for fixture provenance and
[NOTICE.md](NOTICE.md) for source/dependency provenance. Original public-domain
notices and the upstream Unlimited license are preserved in [LICENSE](LICENSE)
and the unchanged R files under `upstream/`.
