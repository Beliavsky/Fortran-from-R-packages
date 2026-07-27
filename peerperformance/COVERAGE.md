# Computational coverage

Upstream package: `PeerPerformance` 2.4.0.

## Exported routines

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `sharpe` | `sharpe` | Implemented |
| `msharpe` | `modified_sharpe` | Implemented |
| `alphaTesting` | `alpha_testing` | Implemented |
| `sharpeTesting` | `sharpe_testing_asymptotic`, `sharpe_testing_bootstrap` | Implemented |
| `msharpeTesting` | `modified_sharpe_testing_asymptotic`, `modified_sharpe_testing_bootstrap` | Implemented |
| `alphaScreening` | `alpha_screening` | Implemented |
| `sharpeScreening` | `sharpe_screening` | Implemented |
| `msharpeScreening` | `modified_sharpe_screening` | Implemented |
| `targetPeerPerformance` | `target_peer_performance` | Implemented |
| `rollScreening` | `roll_screening` | Implemented |
| `exposureHeterogeneity` | `exposure_heterogeneity` | Implemented |

## Supporting algorithms

- complete-case moment and regression calculations
- factor alpha/beta estimation
- classical OLS and automatic Parzen-HAC covariance
- analytical Sharpe and modified-Sharpe gradients
- IID and circular block bootstraps
- VAR(1)-based data-driven block-size selection
- finite-sample-corrected Storey null proportion
- bootstrap mean-squared-error lambda selection
- Ardia-Boudt positive/zero/negative peer-ratio split
- within-group, cross-group, targeted, and rolling workflows

## Represented differently

- R lists and S3 classes are typed Fortran derived types.
- Fund and factor names are represented by integer positions; callers retain labels.
- R `NA` values are IEEE NaNs.
- R parallel-cluster options are omitted; independent focal-fund loops can be
  parallelized by an application if desired.
- Plot methods are replaced by returned numerical arrays.

## Not compiled

- S3 print, summary, confidence-interval, plotting, and data-frame methods
- graphics and graphics-device code
- `parallel` cluster construction and teardown
- bundled `.rdata` loading
- vignette, package-site, and CRAN release infrastructure

These files remain unmodified in `original/PeerPerformance-2.4.0`.
