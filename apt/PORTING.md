# Porting notes

## Functional mapping

| R routine | Fortran routine | Coverage |
|---|---|---|
| `ciTarFit` | `ci_tar_fit` / `ciTarFit` | Complete numerical fit and both F-tests |
| `ciTarLag` | `ci_tar_lag` / `ciTarLag` | Complete lag path and selection |
| `ciTarThd` | `ci_tar_threshold` / `ciTarThd` | Complete trimmed threshold path |
| `ecmSymFit` | `ecm_symmetric_fit` / `ecmSymFit` | Complete two-equation ECM |
| `ecmAsyFit` | `ecm_asymmetric_fit` / `ecmAsyFit` | Complete split/non-split linear, TAR, and MTAR ECM |
| `ecmAsyTest` | `ecm_asymmetry_tests` / `ecmAsyTest` | Complete H1-H4 restriction set |
| `ecmDiag` | `ecm_diagnostics` / `ecmDiag` | Complete statistics; DW p-value differs as noted below |

The R print, summary, and plot methods perform presentation only. Their underlying coefficients, statistics, search paths, and diagnostic values are fields in the Fortran result types.

## Dependency replacement

- `erer::bsLag` is replaced by explicit, index-safe lag construction.
- `erer::lht` and `car::linearHypothesis` are replaced by a general covariance-based linear F-test.
- `stats::lm` is replaced by dense OLS with LAPACK solves.
- `stats::Box.test(type="Ljung")` is reproduced directly.
- `urca` is imported by the R package but is not called by its executable package code.

## Durbin-Watson p-value

`car::durbinWatsonTest` uses simulation by default, so its p-value varies with R's random stream unless a seed is controlled. The Fortran port always reproduces the Durbin-Watson statistic exactly but reports a deterministic, two-sided large-sample p-value based on the lag-one residual correlation. This difference is explicit in the API.

## Time-series representation

R aligns `ts`, `lag`, `window`, and `ts.union` objects using calendar metadata. Fortran callers provide already aligned equal-length arrays. `start_index` in `ci_tar_fit` is an integer observation index and is used internally to reproduce the common-window logic of `ciTarLag(adjust=TRUE)`.

## Coefficient order

For a split asymmetric ECM with lag order `L`, each equation contains:

1. intercept;
2. positive x changes at lags `1:L`;
3. negative x changes at lags `1:L`;
4. positive y changes at lags `1:L`;
5. negative y changes at lags `1:L`;
6. positive error-correction term;
7. negative error-correction term.

This matches the reordered `db` matrix in the R source and is the ordering assumed by H1-H4.

## Numerical behavior

The model matrices must have full column rank. Singular designs return a nonzero status instead of an R condition. IEEE NaNs are used only in low-level regression fields when a statistic cannot be formed; high-level routines report failure through `status`.
