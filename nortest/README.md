# nortest

Modern, pure, free-form Fortran translation of the five computational functions
in R **nortest 1.0-4**, by Juergen Gross and Uwe Ligges. These are tests of the
**composite normality hypothesis**: the sample mean and sample standard
deviation are estimated from the observations. In particular, the
Anderson-Darling calibration here is not a test against a completely specified
CDF.

## Build and run

Place `nortest` beside `rfortran-core` in the
[Fortran-from-R-packages repository](https://github.com/Beliavsky/Fortran-from-R-packages).
The dependency is reused through `../rfortran-core`; its sources are not
included in this package.

```sh
cd nortest
fpm build
fpm test
fpm run --example
```

Only FPM and a modern Fortran compiler are required. There are no system
BLAS/LAPACK dependencies, compiler extensions, shell commands at runtime, or
OS-specific paths in the Fortran. The source is intended to build with FPM and
gfortran on Windows as well as Unix. Validation used Linux, gfortran 14.2.0,
and FPM 0.13.0; a native Windows run has not been performed.

Do **not** enable `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent
flags: missing-value handling and the documented error results require IEEE
NaN/infinity semantics.

## API

```fortran
use nortest, only: dp, normality_result, normality_success, ad_test, pearson_test
real(dp) :: x(8)
type(normality_result) :: answer
x = [-1.3_dp, -0.8_dp, -0.4_dp, -0.1_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.5_dp]
answer = ad_test(x)
if (answer%status == normality_success) print *, answer%statistic, answer%p_value
answer = pearson_test(x, n_classes = 5, adjust = .false.)
```

All five procedures are pure functions taking a rank-one `real(dp)` sample and
returning `normality_result`. `dp` is imported from `r_kinds` and re-exported.
Inputs are never modified. Every NaN is omitted before the sample-size check,
matching R's removal of `NA` and `NaN`. Infinities are rejected.

| Procedure | Accepted nonmissing sample size | Additional arguments |
|---|---:|---|
| `ad_test(x)` | at least 8 | none |
| `cvm_test(x)` | at least 8 | none |
| `lillie_test(x)` | at least 5 | none |
| `pearson_test(x, n_classes, adjust)` | at least 2 | optional integer class count; optional logical adjustment |
| `sf_test(x)` | 5 through 5000 | none |

`pearson_test` defaults to `ceiling(2 * n**0.4)` classes and `adjust=.true.`.
Its degrees of freedom are `n_classes - 3` when adjusted, otherwise
`n_classes - 1`. Zero degrees of freedom follows R's `pchisq` boundary behavior:
the upper-tail probability is one for a zero statistic and zero for a positive
statistic. Negative degrees of freedom and nonpositive class counts are
rejected. No dense class-count array is allocated, so memory usage is O(n),
even with a large class count.

Result fields:

- `statistic`, `p_value`: computed statistic and the upstream approximate p-value
- `n`: number of non-NaN observations, including any rejected infinity
- `status`: one of the constants below; check this before using numerical fields
- `n_classes`, `df`: Pearson metadata, zero for the other tests
- `p_value_at_floor`: true when AD returns its `3.7e-24` lower approximation
  limit or CvM returns `7.37e-10`; it replaces the CvM warning

| Status constant | Value | Meaning |
|---|---:|---|
| `normality_success` | 0 | numerical result available |
| `normality_invalid_size` | 1 | too few/many nonmissing observations |
| `normality_nonfinite` | 2 | an infinity remains after NaN removal |
| `normality_zero_variance` | 3 | constant sample or numerically degenerate spread |
| `normality_invalid_classes` | 4 | nonpositive classes or negative Pearson degrees of freedom |

On failure, `statistic` and `p_value` are quiet NaNs. Size checks precede
finite-value/variance checks, which precede Pearson class validation. Status
replaces R errors/warnings; no I/O or process termination occurs in library
procedures. The result does not recreate R's `htest` class, names, print methods,
or data-expression metadata.

Very large/small finite samples are rescaled before centering to avoid
squaring overflow/underflow. Shapiro-Francia squared correlation is clamped to
[0, 1] to handle floating-point roundoff. Pearson retains the upstream behavior
where a CDF rounded exactly to one maps beyond the last class and is discarded
by R's `tabulate`. Floating-point bin-edge results can depend on rounding.
R's permissive coercion of fractional or vector-valued `n.classes` is omitted:
the Fortran argument is a scalar integer.

The p-value approximations and statistical limitations of the original
package are retained. In particular, its Pearson adjusted/unadjusted p-values
are approximate bounds rather than an exact fitted-normal null distribution;
see `upstream/man/pearson.test.Rd`.

## Translation coverage

Package status: **substantial**. **5 of 5 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R interface compatibility. The same mappings are recorded in
`fpm.toml` and [API_COVERAGE.md](API_COVERAGE.md).

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `ad.test` | `nortest` | `ad_test` | substantial |
| `cvm.test` | `nortest` | `cvm_test` | substantial |
| `lillie.test` | `nortest` | `lillie_test` | substantial |
| `pearson.test` | `nortest` | `pearson_test` | substantial |
| `sf.test` | `nortest` | `sf_test` | substantial |

No exported computational function is untranslated. R object construction and
presentation are intentionally omitted.

## Verification

`fpm test` runs 795 deterministic assertions, including 104 evaluations against
frozen upstream R 4.5.0 references (102 numerical results and two expected
sample-size errors). Fixtures cover minimum sizes, n=100/101 Lilliefors
branches, all four AD/CvM approximation regions plus their lower limits,
normal-shaped/skewed/uniform/tied/bimodal samples, a 2000-observation extreme
outlier, and the Shapiro-Francia upper limit of 5000 observations. All Pearson
options are compared with R. Additional assertions cover missing values,
input immutability, affine invariance, huge/tiny scales, infinity, empty and
constant samples, invalid classes, zero degrees of freedom, and a maximal
integer class count.

Statistic tolerances are `5e-11 * max(1, abs(reference))`; ordinary p-values
use absolute tolerance `5e-10`. Positive reference tails below `1e-6`, including
subnormal Pearson probabilities, require a positive computed result and an
absolute log-ratio error at most `5e-8`. Expected underflowed zeros and reported
approximation floors are checked exactly.

To regenerate the checked-in fixtures using R (optional, not a build dependency):

```sh
Rscript --vanilla test/reference/generate.R
fpm test --profile debug --flag "-fcheck=all -fbacktrace"
fpm test --profile release
fpm run --example --profile debug
fpm run --example --profile release
```

## License and provenance

GPL-2.0-or-later, preserving the upstream `GPL (>= 2)` license. See
[LICENSE](LICENSE), [NOTICE.md](NOTICE.md), and the unmodified package-owned
R sources/documentation under `upstream/`. The shared `rfortran-core` dependency
is separately MIT-licensed and is not vendored here.
