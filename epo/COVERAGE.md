# Computational coverage

## Summary

The upstream package has one numerical implementation, `epo_`, with S3
wrappers for matrices, tibbles, and `xts` objects. The complete numerical
algorithm is translated.

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `epo` | `epo_optimize`, `epo_from_covariance` | Complete numerical interface |
| `epo.matrix` | `epo_optimize` | Complete |
| `epo.tbl` | Caller supplies a numeric Fortran array | R adapter excluded |
| `epo.xts` | Caller supplies a numeric Fortran array | R adapter excluded |
| `epo.default` | Typed validation and `epo_result` diagnostics | Replaced |
| `epo_` | `epo_optimize`, `epo_from_covariance` | Complete |
| `tbl_to_mtx` | Not required | R/tidyverse adapter excluded |

## Numerical functionality

- Sample means and unbiased sample covariance using denominator `n - 1`.
- Covariance-to-correlation conversion.
- Correlation shrinkage toward identity.
- Reconstruction of the shrunk covariance matrix using original variances.
- Simple EPO weights.
- Anchored EPO weights.
- Endogenous scaling coefficient from upstream footnote 13.
- Exogenous risk-aversion scaling.
- Optional full-investment normalization.
- Direct covariance-input calculations.
- Typed intermediate matrices and diagnostics.

## Not applicable to compiled Fortran

- S3 method dispatch.
- Tibble numeric-column selection.
- `xts` conversion.
- `rlang` and `assertthat` error objects.
- R documentation, website, and package infrastructure.

The original source remains available under `original/` for provenance.
