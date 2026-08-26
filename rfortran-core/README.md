# rfortran-core

`rfortran-core` provides small, tested modern Fortran modules for operations
that recur across translations of R packages. It is intentionally independent
of BLAS, LAPACK, operating-system APIs, and package-specific status types.

The initial modules are:

- `r_kinds`: common real and integer kinds and mathematical constants.
- `r_status`: status values returned by shared numerical procedures.
- `r_missing`: elemental IEEE missing/non-finite predicates.
- `r_descriptive`: missing-aware weighted and unweighted means, variances,
  standard deviations, covariances, correlations, and counts.
- `r_distributions`: stable normal density, distribution, and quantile
  functions, including upper-tail and log-probability forms.
- `r_ordering`: stable ordering indices and in-place value sorting.
- `r_sorting`: configurable average ranks, plus compatibility
  exports for the original shared quantile procedures.
- `r_quantiles`: R type-7 quantiles and medians, inverse weighted-ECDF quantiles,
  cumulative-weight interpolation, the fitdistrplus frequency-target rule,
  isotone boundary interpolation, and the twelve survey quantile rules.
- `r_robust`: median absolute deviation with R-compatible defaults.
- `r_stability`: stable log-sum-exp and log-mean-exp reductions.
- `r_special`: positive-domain digamma, trigamma and log-beta functions,
  regularized gamma and beta functions, plus integer log-factorial and
  log-binomial-coefficient helpers.
- `r_rolling`: valid-only and right-aligned trailing sums, minima, maxima,
  means, variances, standard deviations, covariances, and correlations.
- `r_time_series`: biased autocovariance and autocorrelation sequences.
- `r_transforms`: stable `log1p`, `expm1`, `log1mexp`, softplus, logistic, and
  logit transforms.
- `r_vectors`: R-style lagged and repeated differencing for vectors and matrices.
- `r_mod`: a compatibility facade that re-exports the initial API.

New translations should import narrowly from the defining module:

```fortran
use r_kinds, only : dp
use r_descriptive, only : r_mean, r_variance
```

The `r_mod` facade is intended to ease migration of existing translations.

## Numerical policies

- `r_mean`, `r_variance`, and `r_sd` propagate NaNs by default.
- Passing `na_rm=.true.` omits NaNs, corresponding to R's `na.rm = TRUE`.
- Empty reductions and variance with insufficient degrees of freedom return
  an IEEE quiet NaN.
- `r_variance` uses sample variance (`ddof=1`) by default.
- Paired covariance and correlation require conforming vector shapes. They
  propagate non-finite values by default; `na_rm=.true.` omits pairs containing
  NaNs, while `finite_only=.true.` omits pairs containing any non-finite value.
  Covariance uses `ddof=1` by default, and undefined correlation returns a quiet
  NaN, including when either series is constant.
- Weighted statistics require conforming finite, nonnegative weights with a
  positive retained total. Zero-weight observations are excluded before data
  missingness is considered, matching R's `weighted.mean` behavior.
- Weighted variance and covariance divide by total weight by default. Passing
  `unbiased=.true.` uses the reliability-weight correction
  `sum(w) - sum(w**2)/sum(w)`, matching `cov.wt(..., method="unbiased")` and
  reducing to the usual sample denominator for equal positive weights.
- Weighted pairwise statistics apply `na_rm` and `finite_only` jointly to each
  observation pair. Undefined results and invalid weights return a quiet NaN.
- Autocovariances use the R `acf(..., type="covariance")` divisor `n` at every
  lag. Autocorrelations divide these values by the lag-zero autocovariance.
- Time-series routines reject NaNs and infinities and report a status rather
  than terminating the program.
- Matrix time-series operations treat rows as observations and columns as
  series. Differencing is performed along the first dimension.
- `r_quantile_type7` implements R's default quantile interpolation. NaNs are
  propagated unless `na_rm=.true.` is requested.
- `r_median` is the type-7 quantile at probability 0.5. `r_mad` uses that
  median as its default center and R's default scale factor 1.4826.
- `r_weighted_quantile_ecdf` returns the smallest retained value whose
  cumulative normalized weight reaches the requested probability. Its name is
  explicit because the other weighted estimators have different interpolation
  rules.
- `r_weighted_quantile_linear_cdf` linearly interpolates values between adjacent
  normalized cumulative weights, matching the convention migrated from quarks.
- `r_weighted_quantile_frequency_type7` preserves fitdistrplus's
  `1 + (sum(weights) - 1) * probability` target and exact-boundary rule. It is
  not equivalent to expanding integer frequencies and applying base R type 7
  when a target falls inside a repeated-value block.
- `r_weighted_quantile_isotone` preserves isotone's weighted average of the two
  adjacent values when the target exactly equals a cumulative-weight boundary.
- `r_weighted_quantile_survey` implements the math, school, Shah-Vaish, and
  Hyndman-Fan 1--9 rules selected by the public `r_qrule_*` constants.
- All weighted quantile procedures require conforming arrays, a finite
  probability in `[0, 1]`, finite nonnegative weights, and a positive retained
  total. Zero-weight observations are ignored before data missingness is
  considered; invalid inputs return a quiet NaN.
- Average ranks use exact ties by default; callers may supply a relative tie
  tolerance when preserving an existing numerical API requires it.
- Normal lower and upper tails use `erfc`; logarithmic extreme tails use an
  asymptotic expansion to avoid underflow. `r_qnorm` accepts ordinary or log
  probabilities and lower or upper tails, and refines positive values against
  the small upper tail to avoid cancellation.
- Digamma and trigamma currently accept finite positive arguments. Invalid
  arguments return an IEEE quiet NaN; package wrappers may preserve older
  package-specific error conventions.
- Regularized gamma P/Q and beta accept positive finite shape parameters and
  return quiet NaNs for invalid domains. Package wrappers may retain stricter
  validation or historical endpoint sentinels.
- Integer log-factorials require `n >= 0`; integer log-choose requires
  `0 <= k <= n`. Invalid factorials return a quiet NaN and impossible choices
  return negative infinity. Small symmetric `k` values use a stable sum that
  avoids cancellation between large log-gamma values.
- Rolling means use complete trailing windows. `r_roll_mean_valid` returns
  only complete windows; `r_roll_mean_right` returns the input length with
  leading quiet NaNs. A NaN affects only windows that contain it.
- Rolling variance and standard deviation use `ddof=1` by default, matching
  R's sample variance and standard deviation. Callers may select another
  nonnegative `ddof` and may request per-window NaN removal.
- Rolling covariance and correlation apply the paired scalar policies within
  each complete trailing window. The paired inputs must have equal lengths;
  right-aligned results contain leading quiet NaNs.
- Rolling sums, minima, and maxima propagate NaNs by default. Callers may omit
  NaNs or all non-finite values. Empty filtered windows follow R reduction
  identities: zero for sums, positive infinity for minima, and negative
  infinity for maxima. Finite-input rolling sums use a linear-time sliding
  update.
- `r_log1p` follows the real domain `x >= -1`; `r_log1mexp` follows `x <= 0`.
  Values outside those domains return a quiet NaN and boundary values return
  negative infinity. `r_log1pexp` avoids overflow for large positive inputs.
- `r_logistic` is stable in both tails. `r_logit` returns infinities at zero
  and one and a quiet NaN outside the probability interval.
