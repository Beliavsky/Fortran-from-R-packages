# mitools

Modern free-form Fortran translation of the computational code in the R package
`mitools` 2.4 by Thomas Lumley.

The upstream package provides tools for analyzing multiply imputed data sets.
Its main numerical routine is `MIcombine`, which combines estimates and
variance-covariance matrices with Rubin's rules and computes degrees of freedom
and the fraction of missing information. The remaining upstream code is mostly
R data-frame, formula, expression-evaluation, S3, and database infrastructure.

## Implemented numerical API

The public `mitools` module exports:

- `mi_combine`: scalar and multivariate Rubin-rule combination of estimates and
  within-imputation variances/covariances.
- `mi_standard_errors`: standard errors from the combined covariance matrix.
- `mi_confidence_intervals`: Student-t confidence intervals using the combined
  degrees of freedom.
- `mi_summary`: numerical equivalent of the table produced by
  `summary.MIresult`, including optional exponentiation of log effects.
- `imputation_list` and helpers for numeric imputation cubes, extraction,
  dimensions, row binding, and column binding.
- `pv_select` and `pv_materialize` for iterating over one or more sets of
  plausible values without R expression rewriting.

The derived type `mi_result` retains the computational fields of the upstream
`MIresult`: `coefficients`, `variance`, `nimp`, `df`, and `missinfo`.

## Shared dependency

This package uses the repository-local `rfortran-core` package through

```toml
rfortran-core = { path = "../rfortran-core" }
```

and imports `dp` from `r_kinds`, covariance from `r_descriptive`, and the
central Student-t quantile from `r_distributions`. No dependency source is
copied into this package.

## Build

Place this directory at the repository root next to `rfortran-core`, then run:

```text
fpm build
fpm test
fpm run --example mi_example
```

No system BLAS, LAPACK, or ARPACK library is required.

## Array conventions

For `mi_combine`, estimates have shape `(n_parameter, n_imputation)` and
variance-covariance matrices have shape
`(n_parameter, n_parameter, n_imputation)`.

An `imputation_list` stores numeric data as
`(n_row, n_column, n_imputation)`.

Plausible values are stored as `(n_row, n_variable, n_replicate)`.

## Scope

Database connections, DBI/ODBC queries, S3 dispatch, formula parsing,
`model.frame`, expression rewriting/evaluation, and arbitrary R data-frame
updates are intentionally not translated. They are R-specific interfaces rather
than numerical algorithms. See `API_COVERAGE.md` for the detailed mapping.

## License and provenance

The upstream `DESCRIPTION` records `License: GPL-2`. This translation includes
GNU GPL version 2 in `LICENSE`. Upstream metadata and the three R source files
are retained under `reference/upstream/` for provenance. See `NOTICE.md` and
`PROVENANCE.md` for details.
