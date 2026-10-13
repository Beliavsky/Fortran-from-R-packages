# Kendall

Modern, pure, free-form Fortran translation of the computational API of
A.I. McLeod's **Kendall 2.2.2** R package. It provides Kendall tau-b,
Mann-Kendall monotonic-trend tests, and the seasonal Mann-Kendall test.
All maintained real arithmetic uses the shared `dp` kind (IEEE binary64).

## Build and use

Place this directory at the repository root, beside `rfortran-core`:

```text
Fortran-from-R-packages/
  Kendall/
  rfortran-core/
```

From `Kendall`:

```sh
fpm build
fpm test
fpm run --example kendall_demo
```

FPM and gfortran are the only build requirements, including on Windows.
There are no system BLAS/LAPACK links, copied dependencies, R runtime calls,
or operating-system-specific calls in the library or example. Linux/gfortran
was tested; Windows was not independently executed.

```fortran
use kendall_module, only : dp, kendall, kendall_result, kendall_ok
real(dp) :: x(4), y(4)
type(kendall_result) :: ans

x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
y = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp]
call kendall(x, y, ans)
if (ans%status == kendall_ok) print *, ans%tau, ans%sl
```

## API

All three public procedures are pure subroutines in `kendall_module`.
The module re-exports `dp` from `r_kinds`.

- `kendall(x, y, result [, method])`: association of equal-length real vectors
- `mann_kendall(x, result [, method])`: association of observation positions
  with the ordered observations
- `seasonal_mann_kendall(x, period, result)`: aggregate tests for `x(s::period)`
  for `s = 1, ..., period`; `period` is a positive integer

The `kendall_result` fields `tau`, `sl`, `s`, `d`, and `vars` correspond to
R's `tau`, `sl`, `S`, `D`, and `varS`. Additional fields are `n_used`, `status`,
and `exact`. Always inspect `status`; no routine prints, stops, or raises an R
warning. `n_used` counts complete observations (zero on unequal lengths or an
invalid method). Invalid inputs leave statistics NaN. A constant variable has
`s = d = vars = 0`, but undefined `tau` and `sl` remain NaN. An unsupported exact
request preserves computable moments and tau, with `sl = NaN`.

| Status constant | Value | Meaning |
| --- | ---: | --- |
| `kendall_ok` | 0 | Successful test |
| `kendall_size_mismatch` | 1 | Association vectors have different lengths |
| `kendall_too_few` | 2 | Fewer than three complete observations, overall or in a season |
| `kendall_nonfinite` | 3 | A retained observation is infinite |
| `kendall_degenerate` | 4 | Constant variable, or all seasonal contributions are zero |
| `kendall_invalid_method` | 5 | Unrecognized method selector |
| `kendall_exact_ties` | 6 | Exact method requested with ties in either variable |
| `kendall_exact_limit` | 7 | Exact method requested for more than 100 complete pairs |
| `kendall_invalid_period` | 8 | Seasonal period is zero or negative |

### Missing values and seasonal layout

NaN represents a missing observation. `kendall` omits a pair when either member
is NaN; infinities in retained pairs are rejected. Comparisons do not subtract
values, so extreme finite observations cannot overflow while scoring pairs.
Equality of binary64 values determines ties, without tolerances or rounding.
Signed zeros are equal. R factor coercion and S3/time-series classes are omitted.

Mann-Kendall uses original observation positions. Seasonal data must retain
missing time slots as NaN; deleting slots would change seasonal membership.
Each season needs at least three non-NaN values. Unequal season lengths and a
partial final cycle are supported. Constant seasons contribute zero score,
denominator and variance. If all seasons are constant, the aggregate is
undefined. The starting season merely relabels groups and does not affect their
aggregate, so the API needs no `start` argument or R `ts` metadata.

### Statistical methods and deliberate differences

The default `method = kendall_asymptotic` applies the upstream runtime's
continuity-corrected normal approximation, including when there are no ties.
It uses `max(0, abs(S) - 1) / sqrt(varS)`. Seasonal aggregation instead uses
`abs(sum(S)) / sqrt(sum(varS))`, without continuity correction, exactly as the
R seasonal wrapper does.

The upstream documentation promises an exact no-ties test, but its actual
Fortran branch is unreachable. The condition `if (SW)` appears inside
`if (.not. SW)`, and `SW` actually identifies constant variables, despite its
comment. This port does not claim exact-p compatibility with that broken branch.

`method = kendall_exact` is a documented **extension** for no-ties data with
3 through `kendall_exact_max_n = 100` complete observations. It computes
`P(abs(S_random) >= abs(S_observed))` by the inversion-count recurrence over
permutations, with inclusive tails. It uses normalized, nonnegative probability
sums in binary64; "exact" describes the finite-sample distribution, not
arbitrary-precision arithmetic. Ties and sizes above the limit return an
explicit status, without silently switching methods.

Other deliberate corrections:

- Three complete observations receive a computed p-value. Upstream computes
  moments but returns p = 1 and error 12 for every three-point case.
- Constant variables yield an explicit undefined result, rather than upstream's
  misleading tau = 1 and error 12.
- Complete-case counts are validated after omission. Upstream checks the original
  length, then may return placeholder values for fewer than three complete cases.
- Binary64 observations do not collapse to ties through upstream's binary32
  `as.single` conversion. Consequently, very close values can produce different
  scores as well as more accurate statistics.
- Shared `r_pnorm` computes the small upper tail directly. This avoids upstream's
  subtractive cancellation and its hard normal-tail cutoff. For an increasing
  sequence of length 100, upstream returns p = 0; the default here returns
  approximately `3.6371896876271862e-49`.

The score/denominator/tie-corrected variance formulas agree with upstream on
ordinary valid inputs, with binary64 numerical improvements. The seasonal
variance is the sum of within-season variances; neither upstream nor this port
adjusts for serial correlation or cross-season covariance. As with the R
package, correlated observations can invalidate the nominal significance level;
bootstrap workflows belong to callers and are not translated from examples.

Pair scoring is O(n²), with O(n) working memory and O(n log n) tie sorting.
The optional bounded exact recurrence has O(n⁴) worst-case time and O(n²)
memory. Do not use `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent
flags: NaN/Infinity handling is part of the documented API.

## Translation coverage

Package status: **substantial**. **3 of 3 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility; each mapping remains substantial because R
interfaces and the numerical differences above are explicit.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `Kendall` | `kendall_module` | `kendall` | substantial |
| `MannKendall` | `kendall_module` | `mann_kendall` | substantial |
| `SeasonalMannKendall` | `kendall_module` | `seasonal_mann_kendall` | substantial |

There are no untranslated exported computational functions. Printing, summary
presentation, datasets, R classes and package-initialization helpers are outside
the denominator. The imported `boot` package is used by R documentation examples,
not by these three computational functions. No `boot` dependency is needed.
See [API_COVERAGE.md](API_COVERAGE.md) and the machine-readable mapping in
[fpm.toml](fpm.toml).

## Verification

`fpm test` includes book examples, exhaustive no-ties enumeration for n = 3..8,
all 7! permutations for three tied-variance configurations, seasonal sums,
missing slots, constant seasons, short inputs, invalid arguments, precision
preservation, extreme finite values and very small tails.

Optional independent R comparison, from this directory with R, Kendall 2.2.2
and gfortran installed:

```sh
Rscript validation/run_reference.R
```

It compiles maintained sources and the sibling dependency into a temporary
directory. R is not needed for normal builds/tests. The deterministic script
checks 240 association/trend cases, 31 seasonal cases, and 96 exact-distribution
checks against R's `cor.test` or analytical extreme tails. See
[validation/RESULTS.md](validation/RESULTS.md) for tested versions and results.

## License and provenance

GPL-2.0-or-later, preserving the upstream license. [LICENSE](LICENSE) includes
the GPL v2 text; [NOTICE.md](NOTICE.md) identifies upstream authors, exact revision,
source paths, preserved notices and shared dependencies. Archived upstream
sources are outside FPM's source directories; the unchanged fixed-form native
source has suffix `.f.txt` so it cannot be compiled accidentally. All maintained
Fortran, including tests and validation, is free-form.
