# poLCA: modern Fortran translation

This directory translates the computational core of the R package **poLCA 1.6.0.2**, by Drew A. Linzer and Jeffrey Lewis, to modern free-form Fortran with FPM packaging. poLCA estimates latent class models and latent class regression models for polytomous manifest variables using expectation-maximization and Newton-Raphson updates.

The translation is numerical rather than an emulation of R. Callers pass integer manifest-response matrices directly (`1..K_j`, with `0` for missing responses) and, for latent class regression, a numeric predictor matrix that normally includes an intercept column. R formula parsing, data frames, S3 presentation, plotting, and interactive console behavior are intentionally omitted.

## Implemented numerical scope

The public `polca` module exports:

- `polca_fit` for one- or multi-class latent class models, optional latent-class regression, missing manifest responses, deterministic seeded random starts, AIC/BIC, Pearson/deviance goodness-of-fit statistics, posterior modal classes, and poLCA-style outer-product/delta-method standard errors;
- direct translations of the native likelihood, posterior, response-probability, multinomial-prior, and beta-derivative kernels;
- `polca_posterior`, `polca_predcell`, and `polca_entropy`;
- `polca_table_1d` and `polca_table_2d` for expected marginal/conditional tables;
- `polca_simdata` and `polca_rmulti`;
- `polca_reorder_probs`, `polca_coef`, and `polca_vcov`.

The implementation is self-contained. It does **not** link BLAS, LAPACK, MASS, scatterplot3d, `rfortran-compat`, or another translated R package. The symmetric generalized inverse needed by poLCA's score-based standard errors is implemented locally with a Jacobi eigendecomposition.

## Build and test

From the `poLCA` directory:

```text
fpm build
fpm test
fpm run --example latent_class_example
```

No separately installed BLAS or LAPACK is required.

## Minimal example

```fortran
use polca, only : dp, polca_model, polca_fit

type(polca_model) :: fit
integer :: y(8, 3)
integer :: n_choices(3)
real(dp) :: start(2, 2, 3)

! Fill y with one-based categorical responses and start with valid
! class-by-category-by-item probabilities.
n_choices = 2
call polca_fit(y, n_choices, 2, fit, probs_start=start)
```

See `example/latent_class_example.f90` for a complete deterministic example.

## Translation coverage

Package status: **substantial**.

**10 of 10 (100%)** exported computational R functions/S3 computational methods have meaningful Fortran mappings. This fraction measures mapped computational functions, **not** complete compatibility with R formulas, S3 objects, R RNG streams, factor labels, warnings, or presentation behavior.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `poLCA` | `polca_core` | `polca_fit` | substantial |
| `poLCA.simdata` | `polca_api` | `polca_simdata` | substantial |
| `poLCA.reorder` | `polca_api` | `polca_reorder_probs` | complete |
| `poLCA.table` | `polca_api` | `polca_table_1d`, `polca_table_2d` | substantial |
| `poLCA.entropy` | `polca_api` | `polca_entropy` | complete |
| `poLCA.predcell` | `polca_api` | `polca_predcell` | complete |
| `poLCA.posterior` | `polca_api` | `polca_posterior` | substantial |
| `rmulti` | `polca_api` | `polca_rmulti` | substantial |
| `coef.poLCA` | `polca_api` | `polca_coef` | complete |
| `vcov.poLCA` | `polca_api` | `polca_vcov` | complete |

The coverage denominator excludes the presentation-only `print.poLCA` and `plot.poLCA` S3 methods and internal helpers. Exact notes for every mapping are recorded in `fpm.toml` and `API_COVERAGE.md`.

## Validation

The deterministic Fortran suite checks categorical sampling, multinomial priors, missing responses, one-class empirical probabilities, multi-class EM fitting, latent-class regression, standard errors, prediction, entropy, expected tables, class reordering, coefficient/covariance accessors, and simulation.

An independent Python/SciPy parity check optimizes the same observed-data likelihood directly rather than reusing the EM algorithm. At the validated checkpoint it agrees with Fortran to near floating-point precision for both a two-class LCA and a two-class latent class regression, including fitted response probabilities, class shares, regression coefficients, and score/delta-method standard errors. See `VALIDATION.md`.

## Upstream provenance and citation

The upstream package is:

- Drew A. Linzer and Jeffrey Lewis, **poLCA: Polytomous Variable Latent Class Analysis**, R package 1.6.0.2.
- Linzer, Drew A., and Jeffrey Lewis (2011), "poLCA: an R Package for Polytomous Variable Latent Class Analysis," *Journal of Statistical Software* 42(10), 1-29.

The original package is licensed under **GPL version 2 or later**. This translation remains under that license. See `LICENSE`, `NOTICE.md`, and the preserved source material in `upstream/`.
