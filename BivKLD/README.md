# BivKLD

Modern free-form Fortran translation of the computational core of the R package
**BivKLD 0.1.0**, with an FPM package interface.

The upstream package estimates directed bivariate Kullback-Leibler divergence
using Gaussian kernel density estimates and also supplies exact formulas for
several model families. The Fortran translation preserves those numerical
operations while omitting plotting, R data-frame/list dispatch, S3 printing,
attributes, and other R-specific presentation behavior.

## Build

Place this directory at the root of `Fortran-from-R-packages` (or build it by
itself; this translation has no external FPM dependencies), then run:

```text
fpm build
fpm test
fpm run --example basic
fpm run --example pairwise
```

No system BLAS or LAPACK library is required. The implementation uses explicit
closed-form two-by-two linear algebra and direct Gaussian kernel evaluation.
Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or similar
options if you need the documented NaN/infinity behavior.

## Public API

Use the facade module:

```fortran
use bivkld
```

The translated R operations are:

- `biv_kld`: directed two-sample kernel divergence. The Fortran procedure is a
  subroutine with scalar `estimate` output and optional `info` status. Optional
  `bivkld_estimate` output carries bandwidths, grid axes, densities, cell area,
  standardization, and selector information.
- `biv_kld_matrix`: all directed pairwise divergences for an array of
  `biv_sample` values on one common grid.
- `biv_kld_discrete`: generic interface for vector and rank-2 probability
  tables; higher-rank R arrays are not separately overloaded.
- `biv_kld_normal`, `biv_kld_pareto2`, and
  `biv_kld_independent_weibull`: analytic formulas.
- `hns_bandwidth`, `hscv_bandwidth`, and `kde_density`: public lower-level
  helpers useful for numerical workflows.

Kernel routines report these status codes: `BIVKLD_SUCCESS = 0`,
`BIVKLD_INVALID_INPUT = 1`, `BIVKLD_INVALID_BANDWIDTH = 2`,
`BIVKLD_NUMERICAL_FAILURE = 3`, and `BIVKLD_INVALID_OPTION = 4`. Exact
functions return IEEE quiet NaN for input-validation failures, matching R's
error cases without terminating a Fortran program; discrete support mismatch
returns positive infinity as in R.

### Bandwidth selectors

`bandwidth="normal"` implements the bivariate normal-scale `Hns` rule.
`bandwidth="scv"` minimizes the unbinned bivariate smoothed-cross-validation
criterion with a full symmetric positive-definite bandwidth matrix. The SCV
criterion follows the Gaussian pair-functional form used by `ks::Hscv`, but
this translation uses a streamlined normal-scale pilot and a three-parameter
Nelder-Mead search rather than reproducing all `ks` pilot, pre-scaling,
binning, and optimizer options. This is a material compatibility difference
and is why the kernel mappings are marked `substantial` rather than `complete`.

The repository's existing `ks` translation was evaluated before implementing
this code. It does not currently translate multivariate `Hscv`, and its elected
GPL-2.0-only package license is not used as a linked dependency of this GPL-3
translation. No `ks` source is copied or vendored here; the small KDE and
bandwidth operations needed specifically by BivKLD are independently
implemented from the mathematical definitions and upstream behavior.

## Translation coverage

Package status: **substantial**.

Coverage: **6 of 6 (100%)** exported computational R functions have meaningful
Fortran mappings. This fraction measures mapped computational functions, not
complete R interface or behavioral compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `biv_kld` | `bivkld_kernel` | `biv_kld` | substantial |
| `biv_kld_matrix` | `bivkld_kernel` | `biv_kld_matrix` | substantial |
| `biv_kld_discrete` | `bivkld_exact` | `biv_kld_discrete_vector`, `biv_kld_discrete_matrix` | substantial |
| `biv_kld_normal` | `bivkld_exact` | `biv_kld_normal` | complete |
| `biv_kld_pareto2` | `bivkld_exact` | `biv_kld_pareto2` | complete |
| `biv_kld_independent_weibull` | `bivkld_exact` | `biv_kld_independent_weibull` | complete |

See `API_COVERAGE.md` for compatibility notes and `NOTICE.md` for upstream
provenance and licensing information.
