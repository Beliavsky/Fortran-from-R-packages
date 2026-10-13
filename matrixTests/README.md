# matrixTests: modern Fortran

Numerical translation of CRAN matrixTests **0.2.3.1**, by Karolis Koncevičius.
This top-level FPM package supplies all 21 row/column test families: t tests,
variance-ratio F tests, Pearson correlation, ordinary/Welch one-way ANOVA,
Kruskal-Wallis, van der Waerden, Bartlett, Fligner-Killeen, Levene,
Brown-Forsythe, Wilcoxon tests, Jarque-Bera, Anderson-Darling, cosinor, and
exact/asymptotic two-sample Kolmogorov-Smirnov.

## Build and run

Place this directory as `matrixTests/` at the root of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages).
The sibling packages `rfortran-core`, `matrixStats`, `exactRankTests`, `nortest`,
and `moments` must be present. They are reused by FPM path dependencies and are
not copied here.

```sh
fpm build
fpm test
fpm run --example --all
```

A tested, revision-pinned dependency on the
[Beliavsky fork of fortran-lapack](https://github.com/Beliavsky/fortran-lapack)
provides the pure least-squares call chain. It supplies its own Fortran BLAS and
LAPACK. No separately installed BLAS/LAPACK, platform-specific shell calls,
or system library link directives are required. FPM can fetch it on Windows
with gfortran and Git. This package was verified on Linux/gfortran; Windows
execution was not available. Never compile with `-ffast-math`, `-Ofast`, or
`-ffinite-math-only`: missing-value and infinity behavior requires IEEE arithmetic.

FPM's version parser accepts three numeric components, so the manifest's
translation version is `0.2.3`; the authoritative R source version is `0.2.3.1`.

## API

Use module `matrixtests`, which exports `dp`, `test_result`, status constants,
and all 42 procedures listed below. Every public row/column function returns
an allocatable array of `test_result`, one element per tested row or column.
All real input arrays and numerical results use the single `dp` kind imported
from `r_kinds`.

```fortran
use matrixtests, only : dp, test_result, row_t_welch
real(dp) :: x(2, 6), y(2, 4)
type(test_result), allocatable :: results(:)
! Initialize x and y, then:
results = row_t_welch(x, y, alternative='greater', conf_level=0.95_dp)
```

- Matrix observations run along columns for `row_*` and rows for `col_*`.
- Two-sample `y` may have one tested row/column, broadcast over `x`; otherwise
  its number of tests must equal that of `x`. Paired and correlation tests
  require equal numbers of observations.
- Group tests take an integer vector `g` of observation labels. Positive values
  are arbitrary groups; zero/negative labels omit those observations. Convert
  R factors/strings to integers before calling.
- Options are scalar and common to all tests in a call. `null` defaults to zero
  (one for `f_var`) and must be finite; a variance null must be positive.
  `alternative` accepts exactly `two.sided`, `less`, or `greater`.
  `conf_level` defaults to 0.95, accepts [0,1], and IEEE NaN disables intervals.
  Use several calls to supply different options to different tests.
- `exact` is an optional logical. If absent, Wilcoxon uses sample sizes below
  50 and KS uses the product of sample sizes below 10,000. Wilcoxon zeros and
  ties disable exact computation as in upstream. `correct` defaults to true
  for approximate Wilcoxon inference. The returned `exact`/`corrected` flags
  report the method selected (even when an empty sample makes inference undefined).
- Cosinor takes a time vector `t` and optional positive finite `period` (24 by
  default). Acrophase is reported in time units in [0,period).

The result type is a union of the numerical columns for all families.
Dots in R column names become underscores: `mean.diff` becomes `mean_diff`,
`conf.low` becomes `conf_low`. The exceptions are `df.between`/`df.model`,
which use `df_num`, and `df.within`/`df.residual`, which use `df_denom`.
All three kinds of null estimates use `null`. Read fields corresponding to
the selected family's upstream columns; additional summary fields may also be
populated. Unassigned fields start as IEEE NaN, and unassigned counts as zero. The inference
status is `mt_ok=0`, `mt_invalid=1`, `mt_insufficient=2`, `mt_degenerate=3`,
`mt_limit=4`, or `mt_nonfinite=5`. Undefined statistics/p-values are NaN.
Diagnostic summaries that remain meaningful can still be populated after a
nonzero status. No warnings, exceptions, or process exits are emitted by the
library for invalid inputs.

Every maintained numerical procedure is pure except the Jarque-Bera adapters:
the existing compatible `moments` API is impure. The test/example programs
perform I/O. All dummy arguments have explicit intent and FORD documentation.

## Missing values and numerical conventions

IEEE NaN represents both R NA and NaN. Each test removes missing observations
independently; paired tests remove missing differences, correlation removes
pairs whose sum is NaN. The diagnostic `obs_x` and `obs_y` in paired tests
still count each original sample independently. Signed-rank counts exclude
zero differences. Bartlett and Welch ANOVA drop groups with fewer than two
observations. R data frames, names, formula methods, recycling, warning text,
character/factor conversion, and vector-valued options are not replicated.

Variance-ratio, Levene/Brown-Forsythe, and cosinor tests omit infinite responses
as upstream. Rank tests and KS retain infinities as ordered values. Other
moment-based tests return nonfinite/undefined diagnostics for infinite input
instead of relying on incidental R missing-value propagation. Overflowing
moment calculations are reported, not claimed as successful inference.

Intentional numerical differences and limits:

1. Direct upper distribution tails preserve small probabilities where upstream
   uses `1-CDF` and rounds to zero. Welch t degrees of freedom use normalized
   variance fractions to avoid squaring large variances. Pearson, Jarque-Bera,
   and the reused Anderson-Darling implementation use shift/scale safeguards.
2. Exact Wilcoxon reuses `exactRankTests` integer dynamic programming, limited
   by its default 5,000,000-cell workspace. An over-limit exact request returns
   `mt_limit`; it does not silently use a normal approximation.
3. Exact KS uses conditional random-label probabilities at unique pooled-value
   boundaries, including ties, with an explicit `n_x*n_y <= 25,000,000` work
   limit. It evaluates upper tails directly. Forced exact requests above the
   limit return `mt_limit`. Its two-sided asymptotic theta series is converged
   more accurately than R's truncation: for n_x=n_y=50, D=.18, the result is
   .39273070794065437 versus R 4.5.0's .39273380987443274. Repeated infinities
   are treated as ties correctly: x=[-Inf,-Inf,0,2], y=[-Inf,1,Inf] gives exact
   two-sided 29/35 (independent label enumeration), where R 4.5.0 gives 31/35.
4. Cosinor uses pure pivoted LAPACK minimum-norm least squares at relative
   rank tolerance 1e-7. Full-rank coefficients and inference match upstream
   within numerical tolerance. Rank-deficient coefficients can differ from
   R's ordered QR, so both cosinor mappings are explicitly partial. Parameters
   may still be returned when inferential degrees of freedom are insufficient.
5. Levene/Brown-Forsythe require at least one group with three observations in
   each individual test. This repairs upstream's matrix-global `all` check,
   which can make one row's eligibility depend on another row.
6. Welch ANOVA with zero-variance groups returns explicit degeneracy rather
   than upstream's unstable arithmetic. Finite null hypotheses are required.

## Validation

`fpm test` runs 1,920 committed actual-R numerical field fixtures, checks every
column adapter against its row counterpart, compares exact KS with 3,270
small tied/untied random-label allocations, and exercises independent
identities, empty/NaN data, zero variances, singleton groups, broadcasting,
confidence endpoints, extreme tails, and very large scales. KS approximate
R fixtures allow the documented upstream theta truncation difference; exact
and independent mathematical checks use strict tolerances.

`tools/generate_reference.R` regenerates `test/reference.txt` using the actual
matrixTests 0.2.3.1 package; `test/REFERENCE_VERSION.txt` records R and matrixStats
versions. R is needed only to regenerate reference records, not to build or test
this Fortran package. Adapter/coverage generation scripts are reproducible
maintenance conveniences; generated Fortran sources are already included.

The release checks are default build/test/examples, checked-debug
build/test/examples with `-fcheck=all`, and release build/test/examples.
See `NOTICE.md` for source pins, authorship, dependencies, and license details.

## Translation coverage

Package status: **substantial**. **42 of 42 (100%)** exported computational R functions have numerical mappings. The fraction measures mapped computational functions, not complete R compatibility. Forty mappings are substantial; the two cosinor mappings are partial because rank-deficient coefficient selection differs from R.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `row_t_equalvar` | `matrixtests` | `row_t_equalvar` | substantial |
| `col_t_equalvar` | `matrixtests` | `col_t_equalvar` | substantial |
| `row_t_paired` | `matrixtests` | `row_t_paired` | substantial |
| `col_t_paired` | `matrixtests` | `col_t_paired` | substantial |
| `row_t_onesample` | `matrixtests` | `row_t_onesample` | substantial |
| `col_t_onesample` | `matrixtests` | `col_t_onesample` | substantial |
| `row_t_welch` | `matrixtests` | `row_t_welch` | substantial |
| `col_t_welch` | `matrixtests` | `col_t_welch` | substantial |
| `row_cor_pearson` | `matrixtests` | `row_cor_pearson` | substantial |
| `col_cor_pearson` | `matrixtests` | `col_cor_pearson` | substantial |
| `row_f_var` | `matrixtests` | `row_f_var` | substantial |
| `col_f_var` | `matrixtests` | `col_f_var` | substantial |
| `row_oneway_equalvar` | `matrixtests` | `row_oneway_equalvar` | substantial |
| `col_oneway_equalvar` | `matrixtests` | `col_oneway_equalvar` | substantial |
| `row_oneway_welch` | `matrixtests` | `row_oneway_welch` | substantial |
| `col_oneway_welch` | `matrixtests` | `col_oneway_welch` | substantial |
| `row_kruskalwallis` | `matrixtests` | `row_kruskalwallis` | substantial |
| `col_kruskalwallis` | `matrixtests` | `col_kruskalwallis` | substantial |
| `row_waerden` | `matrixtests` | `row_waerden` | substantial |
| `col_waerden` | `matrixtests` | `col_waerden` | substantial |
| `row_bartlett` | `matrixtests` | `row_bartlett` | substantial |
| `col_bartlett` | `matrixtests` | `col_bartlett` | substantial |
| `row_jarquebera` | `matrixtests` | `row_jarquebera` | substantial |
| `col_jarquebera` | `matrixtests` | `col_jarquebera` | substantial |
| `row_flignerkilleen` | `matrixtests` | `row_flignerkilleen` | substantial |
| `col_flignerkilleen` | `matrixtests` | `col_flignerkilleen` | substantial |
| `row_levene` | `matrixtests` | `row_levene` | substantial |
| `col_levene` | `matrixtests` | `col_levene` | substantial |
| `row_brownforsythe` | `matrixtests` | `row_brownforsythe` | substantial |
| `col_brownforsythe` | `matrixtests` | `col_brownforsythe` | substantial |
| `row_wilcoxon_onesample` | `matrixtests` | `row_wilcoxon_onesample` | substantial |
| `col_wilcoxon_onesample` | `matrixtests` | `col_wilcoxon_onesample` | substantial |
| `row_wilcoxon_paired` | `matrixtests` | `row_wilcoxon_paired` | substantial |
| `col_wilcoxon_paired` | `matrixtests` | `col_wilcoxon_paired` | substantial |
| `row_wilcoxon_twosample` | `matrixtests` | `row_wilcoxon_twosample` | substantial |
| `col_wilcoxon_twosample` | `matrixtests` | `col_wilcoxon_twosample` | substantial |
| `row_cosinor` | `matrixtests` | `row_cosinor` | partial |
| `col_cosinor` | `matrixtests` | `col_cosinor` | partial |
| `row_andersondarling` | `matrixtests` | `row_andersondarling` | substantial |
| `col_andersondarling` | `matrixtests` | `col_andersondarling` | substantial |
| `row_kolmogorovsmirnov_twosample` | `matrixtests` | `row_kolmogorovsmirnov_twosample` | substantial |
| `col_kolmogorovsmirnov_twosample` | `matrixtests` | `col_kolmogorovsmirnov_twosample` | substantial |
