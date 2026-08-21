# Upstream source map

This file maps the principal upstream `locfit` sources to the v0.1.0 Fortran
translation.  `upstream/locfit-R/` contains the exact source package used.

| Upstream source | Fortran module | v0.1.0 coverage |
|---|---|---|
| `src/lfcons.h` | `locfit_constants` | Constants, kernels, links, families, status codes and residual types |
| `src/math.c`, `src/prob.c` | `locfit_math` | Stable exp/logistic, normal tail/CDF, incomplete gamma and beta support |
| `src/weight.c` | `locfit_kernels` | Kernel catalogue, distances, weights, derivatives, Taylor coefficients, moments and convolutions used by selectors |
| `src/lf_fitfun.c` | `locfit_basis` | Standard polynomial/angular basis and derivative basis |
| `src/family.c` | `locfit_families` | Standard family/link likelihood equations including censored branches |
| `src/lf_nbhd.c` | `locfit_core` | General distance calculation, nearest-neighbor bandwidth and neighborhood weights; not the upstream ordered-1D acceleration structure |
| `src/locfit.c` | `locfit_core` | Core weighted local likelihood Newton solve and Gaussian fast path |
| `src/lf_vari.c` | `locfit_core` | Sandwich covariance construction for observation-sum families |
| `src/fitted.c` | `locfit_diagnostics` | Residual conversion and studentization formulas |
| `R/locfit.r` (`gcv`, `aic`, `cp`) | `locfit_diagnostics` | Criterion formulas |
| `R/locfit.r` (`km.mrl`) | `locfit_diagnostics` | Kaplan-Meier mean residual life |
| `R/locfit.r` (`locfit.robust`, `locfit.quasi`) | `locfit_robust` | Iterative high-level weighting loops |
| `src/band.c` (`compsda`, `widthsj`, `kdecri`, `esolve`, `kdeselect`) | `locfit_bandwidth` | KDE bandwidth criteria and selectors |
| `src/density.c`, `src/dens_odi.c` | `locfit_density` | 1-D KDE and 1-D local density likelihood; uses numerical integration for the local density normalizer |
| `src/ev_interp.c` | `locfit_interpolation` | Linear/Hermite/rectangular-cell interpolation primitives |
| `src/m_*.c` | `locfit_linalg` plus Fortran intrinsics | Dense linear solve, inverse, rank helpers required by translated core; upstream internal matrix ABI is not reproduced |

The R formula parser, S3 object machinery, lattice/base plotting and printing
methods are intentionally outside the scope of the Fortran port.
