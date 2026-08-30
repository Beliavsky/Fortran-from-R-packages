# Upstream computational mapping

This document maps the main numerical source areas in R `vars` 1.6-1 to the
maintained Fortran implementation.

| Upstream R source | Fortran implementation | Notes |
| --- | --- | --- |
| `VAR.R` | `src/vars_regression.f90` | VAR design matrix and least-squares estimation; constant/trend/both/none, centered seasonal dummies, exogenous regressors. |
| `VARselect.R` | `src/vars_regression.f90` | Fixed-sample AIC, HQ, SC/BIC, and FPE lag selection. |
| `restrict.R` | `src/vars_regression.f90` | Manual zero restrictions and sequential elimination by absolute t statistic. |
| `logLik.varest.R` | `src/vars_regression.f90` | Gaussian VAR likelihood. |
| `A*.R`, `B*.R` | `var_model%a`, `var_model%coef`; `structural_impact` | Typed arrays replace R accessor methods. |
| `Phi*.R`, `Psi*.R` | `src/vars_dynamics.f90` | MA and orthogonalized MA recursions. |
| `roots.R` | `src/vars_dynamics.f90` | Companion matrix and general real eigensystem. |
| `predict.*.R`, `.fecov*` in `internal.R` | `src/vars_dynamics.f90` | Forecast recursion and covariance accumulation. |
| `irf*.R`, `fevd*.R`, `.irf` | `src/vars_dynamics.f90` | Reduced and structural IRF/FEVD arrays. |
| `BQ.R` | `src/vars_dynamics.f90` | Blanchard-Quah long-run identification. |
| `vec2var.R` | `vec2var_coefficients` | Numerical VECM-to-level-VAR lag conversion; R `ca.jo` extraction/object construction omitted. |
| `arch.R` | `src/vars_diagnostics.f90` | Univariate and multivariate ARCH-LM statistics. |
| `normality.R` | `src/vars_diagnostics.f90` | Univariate and multivariate Jarque-Bera decomposition. |
| `serial.R` | `src/vars_diagnostics.f90` | Portmanteau, adjusted Portmanteau, Breusch-Godfrey LM and Edgerton-Shukur statistics. |
| `causality.R` | `src/vars_causality.f90` | Homoskedastic Granger Wald/F and instantaneous-causality chi-square calculations. |
| `SVAR.R` | `src/vars_structural.f90` | A/B structural likelihood, scoring iteration, covariance approximation, LR statistic, sign normalization. |
| `SVEC.R` | `src/vars_structural.f90` | Long-run multiplier and structural scoring from numeric Johansen/VECM arrays. |
| `.boot` in `internal.R` | `src/vars_bootstrap.f90` | Residual resampling recursion and reduced-form IRF bootstrap; resampling indices are supplied by caller instead of invoking R RNG. |
| `stability.R` | external translated `strucchange` | Upstream computation is delegated to `strucchange::efp`; not duplicated. |
| `toMlm.R`, coefficient/fitted/residual methods | typed model fields | R model-object adaptation only. |
| plot/print/summary/fanchart methods | omitted | R-specific presentation and interaction. |

## Deliberate interface changes

R accepts formulas, names, `varest`, `vec2var`, `svarest`, `svecest`, and `ca.jo`
objects. Fortran accepts numeric arrays and derived types. This preserves the
computational algorithms while avoiding an R runtime or a clone of R's object
system.

The upstream direct SVAR branch delegates optimization to R's generic `optim`.
The translated package provides the exact structural objective through
`svar_negloglik` and the package-owned scoring estimator through
`svar_fit_scoring`; callers that want another optimizer can minimize the exposed
objective with their preferred Fortran optimization package.

Similarly, robust sandwich-covariance callbacks accepted by R causality methods
are interface/plugin behavior rather than a single `vars` numerical algorithm.
The native causality routines implement the package's homoskedastic core
statistics directly.

## Interface normalization notes

### Structural covariance convention

The R `svarest` presentation object scales its stored residual covariance by 100
and compensates for that convention in downstream R methods.  The Fortran types
store covariance matrices in their natural numerical scale throughout; no display
scaling is applied.

### Granger causality

The native `granger_causality` routine implements the homoskedastic Wald/F test for
an unrestricted fitted VAR.  R formula/model-object overlays for already restricted
equation objects and the optional sandwich `vcov.` callback are interface-level
features and are not cloned.
