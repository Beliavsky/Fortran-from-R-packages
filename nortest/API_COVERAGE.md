# API coverage

Status: **substantial**. **5 of 5 (100%)** distinct exported computational R
functions mapped, with no untranslated functions. The denominator is the exact
five-function `upstream/NAMESPACE` export list. This is a computational mapping
fraction, not a claim of complete R compatibility.

| R function | R source | Public Fortran procedure | Status |
|---|---|---|---|
| `ad.test` | `upstream/R/ad.test.R` | `nortest::ad_test` | substantial |
| `cvm.test` | `upstream/R/cvm.test.R` | `nortest::cvm_test` | substantial |
| `lillie.test` | `upstream/R/lillie.test.R` | `nortest::lillie_test` | substantial |
| `pearson.test` | `upstream/R/pearson.test.R` | `nortest::pearson_test` | substantial |
| `sf.test` | `upstream/R/sf.test.R` | `nortest::sf_test` | substantial |

All public procedures above are implemented in `src/nortest.f90`. The
statistics, fitted-parameter standardization, finite-sample adjustments, and
piecewise p-value calibrations are translated. Missing observations are
removed before checking the required sample size.

Material interface/numerical differences:

- `normality_result` replaces `htest`, named statistics, method strings,
  expression capture, and printing.
- Status codes and NaN numerical fields replace errors on invalid data.
  Infinities, constant samples, and Pearson negative degrees of freedom receive
  explicit statuses instead of R's incidental errors, warnings, or NaNs.
- CvM's lower-accuracy warning becomes `p_value_at_floor`; AD's fixed lower
  approximation is also marked.
- Pearson accepts a scalar integer class count only; no R numeric coercion or
  vector recycling is reproduced. Its zero-degree-of-freedom and CDF-rounding
  boundary behavior are preserved.
- Scaling protects very large and tiny finite values, and the Shapiro-Francia
  squared correlation is clamped to the mathematically valid [0, 1] interval.
- Existing `rfortran-core` implementations provide mean, sample standard
  deviation, correlation, sorting, and normal/chi-square distributions. The
  original R `stats` implementation is not copied.

Original scientific references and test limitations are preserved in
`upstream/man/`. README and machine-readable `fpm.toml` coverage agree.
