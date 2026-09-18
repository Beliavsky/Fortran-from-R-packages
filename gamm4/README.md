# gamm4 Fortran translation

Modern free-form Fortran translation of the computational core of the R package
[`gamm4`](https://cran.r-project.org/package=gamm4), version 0.3-0.

`gamm4` is primarily a bridge between **mgcv** smooth construction and **lme4**
mixed-model fitting. This translation therefore reuses the existing sibling
Fortran translations instead of copying their code.

## Build

Place this directory at the root of `Fortran-from-R-packages` beside `mgcv`,
`lme4`, and `rfortran-linalg`, then run:

```text
fpm build
fpm test
fpm run --example example_gamm4
```

The package has no system BLAS/LAPACK link directives and vendors no dependency
source. `rfortran-linalg` manages its own FPM-compatible linear-algebra
implementation/dependency chain.

## Numerical API

The public umbrella module is `gamm4`.

```fortran
use gamm4
```

Construct each smooth as a `gamm4_smooth_t` containing an explicit basis matrix
and an mgcv `smooth_spec_t` with one quadratic penalty. Optional ordinary random
effects use lme4 `random_term_t` values. Then call:

```fortran
call gamm4_fit(y, x_fixed, smooths, result, status, &
   random_terms=terms, family=gamm4_family_gaussian)
```

For Gaussian models the default is REML; pass `reml=.false.` for ML. The
translated lme4 family constants select binomial, Poisson, Gamma,
inverse-Gaussian, or negative-binomial fits.

The result contains original-basis GAM coefficients and covariance, per-
coefficient EDF contributions, smoothing parameters, smooth random-effect
standard deviations, GAM-only fitted values, full conditional fitted values,
variance parameters, and mixed-model random-effect modes.

## Translation coverage

Package status: **substantial**.

**1 of 1 (100%)** exported computational R functions is mapped. The fraction
measures mapped computational functions rather than complete R compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `gamm4` | `gamm4_fit_mod` | `gamm4_fit` | substantial |

Formula/model-frame handling, S3 object construction, multi-penalty smooths,
and the optional reticulate/CHOLMOD acceleration route are intentionally not
claimed as translated. See `API_COVERAGE.md` and `PORTING_NOTES.md`.

## Validation

The deterministic tests exercise:

- direct `getVb` marginal covariance calculation;
- Gaussian GAMM fitting with both a smooth and an ordinary grouped random
  intercept;
- Poisson smooth fitting through the PIRLS/Laplace path; and
- explicit rejection of unsupported multi-penalty smooths.

An independent SciPy/NumPy check reconstructs the Gaussian profiled mixed-model
objective and the upstream `getVb` matrix formula. See `VALIDATION.md`.

## License and provenance

GPL-2.0-or-later. See `COPYING`, `LICENSE.md`, `NOTICE.md`, and `upstream/`.
